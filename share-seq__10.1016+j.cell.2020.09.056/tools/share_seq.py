"""SHARE-seq joint RNA/accessibility combinatorial indexing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import MolecularState, Scene, Segment, Workflow, feature
from multimodal_spatial import HybridizationLigationRound, modality_split

TITLE="SHARE-seq"
NOTES="01_share-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma et al., <i>Cell</i> (2020)</a>.'
SUMMARY="Fixed permeabilized cells are first tagmented and reverse-transcribed. Three split-pool hybridization rounds assemble the same barcode chain on ATAC fragments and biotinylated cDNA; blocking oligos sequester free adapters, then one ligation seals the chain before the modalities are separated."
CAVEAT="The constant arms, linkers and ordering modifications come from Supplementary Table S1. The 96 well-specific 8-nt barcode values in each round remain placeholders."
BC1=feature("cell_barcode_round1","cell_barcode","combinatorial",group="cell_barcode",part="round 1")
BC2=feature("cell_barcode_round2","cell_barcode","combinatorial",group="cell_barcode",part="round 2")
BC3=feature("cell_barcode_round3","cell_barcode","combinatorial",group="cell_barcode",part="round 3")
UMI_FEATURE=feature("umi","umi","random")

# Ma et al. Supplementary Table S1 / the authors' SHARE-seq protocol.  Each
# barcode value varies by well; the constant arms, splints and modifications do
# not.  The RT primer and the ligatable Tn5 end deliberately expose the same
# 15-nt arm, allowing the same round-1 plate to label both modalities.
SHARED_ANCHOR = "GTCTCGTGGGCTCGG"
ME = "AGATGTGTATAAGAGACAG"
RT_HANDLE = SHARED_ANCHOR + ME
RT_DT = "T" * 14 + "VN"
TSO_WRITTEN = "AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G"
RNA_PCR = "AAGCAGTGGTATCAACGCAGAGT"
READ1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
BLOCKED_ME_COMP_WRITTEN = "/5Phos/C*T*G* T*C*T* C*T*T* A*T*A* C*A*/3ddC/"

ROUNDS = (
    HybridizationLigationRound(
        "round 1", SHARED_ANCHOR, "CGCGCTGCATACTTG", "CCCATGATCGTCCGA",
        "CCGAGCCCACGAGACTCGGACGATCATGGG", BC1, phosphorylated_5=True),
    HybridizationLigationRound(
        "round 2", "CGCGCTGCATACTTG", "CATCGGCGTACGACT", "ATCCACGTGCTTGAG",
        "CAAGTATGCAGCGCGCTCAAGCACGTGGAT", BC2, phosphorylated_5=True),
    HybridizationLigationRound(
        "round 3", "CATCGGCGTACGACT", "CAAGCAGAAGACGGCATACGAGAT",
        "GTGGCCGATGTTTCG", "AGTCGTACGCCGATGCGAAACATCGGCCAC", BC3),
)

RNA,RNA_P=truseq_library([seg("barcode round 3","C"*8,"cbc",placeholder=True,feature=BC3),seg("barcode round 2","B"*8,"cbc",placeholder=True,feature=BC2),seg("barcode round 1","A"*8,"cbc",placeholder=True,feature=BC1),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI_FEATURE),seg("cDNA","X"*34,placeholder=True)],"SHARE-seq RNA library")
ATAC,ATAC_P=nextera_library([seg("barcode round 3","C"*8,"cbc",placeholder=True,feature=BC3),seg("barcode round 2","B"*8,"cbc",placeholder=True,feature=BC2),seg("barcode round 1","A"*8,"cbc",placeholder=True,feature=BC1),seg("accessible DNA","X"*34,placeholder=True)],"SHARE-seq ATAC library")
FINAL_LIBRARIES=(
    ("Final RNA library",RNA,RNA_P,"Three ligated well barcodes and the RT UMI identify each transcript and cell.","The RNA library's declared primer sites are computed on the completed duplex."),
    ("Final ATAC library",ATAC,ATAC_P,"The same three well barcodes identify the cell of origin of each accessible-DNA fragment.","The ATAC library's declared primer sites are computed on the completed duplex."),)

def oligos():
    return (
        f"<b>Tn5 Read 1 oligo</b> — 5′-{READ1}-3′",
        f"<b>Phosphorylated Read 2 oligo</b> — 5′-/5Phos/{RT_HANDLE}-3′",
        f"<b>Blocked mosaic-end complement</b> — 5′-{BLOCKED_ME_COMP_WRITTEN}-3′",
        f"<b>RT primer</b> — 5′-/5Phos/{RT_HANDLE}[10-nt UMI]/iBiodT/{RT_DT}-3′",
        f"<b>Round 1 barcode</b> — 5′-/5Phos/{ROUNDS[0].barcode_5_arm}[8-nt barcode]{ROUNDS[0].barcode_3_arm}-3′",
        f"<b>Round 1 linker</b> — 5′-{ROUNDS[0].linker}-3′",
        f"<b>Round 1 blocking oligo</b> — 5′-{ROUNDS[0].blocking_oligo}-3′",
        f"<b>Round 2 barcode</b> — 5′-/5Phos/{ROUNDS[1].barcode_5_arm}[8-nt barcode]{ROUNDS[1].barcode_3_arm}-3′",
        f"<b>Round 2 linker</b> — 5′-{ROUNDS[1].linker}-3′",
        f"<b>Round 2 blocking oligo</b> — 5′-{ROUNDS[1].blocking_oligo}-3′",
        f"<b>Round 3 barcode</b> — 5′-{ROUNDS[2].barcode_5_arm}[8-nt barcode]{ROUNDS[2].barcode_3_arm}-3′",
        f"<b>Round 3 linker</b> — 5′-{ROUNDS[2].linker}-3′",
        f"<b>Round 3 blocking oligo</b> — 5′-{ROUNDS[2].blocking_oligo}-3′",
        f"<b>Template-switching oligo</b> — 5′-{TSO_WRITTEN}-3′",
        f"<b>RNA PCR primer</b> — 5′-{RNA_PCR}-3′",
    )


def _initial_rows():
    chromatin = Scene.duplex(
        [Segment("fixed chromatin", "X" * 36, placeholder=True)],
        label="fixed nuclear chromatin")
    mrna = Scene()
    mrna.strand("mRNA", [Segment("transcript", "X" * 24, placeholder=True),
                         Segment("poly(A)", "A" * 18)], label="poly(A) RNA")
    return [*chromatin.rows(), *mrna.rows()]


def _prefix(rounds: int) -> list[Segment]:
    out: list[Segment] = []
    for round_ in reversed(ROUNDS[:rounds]):
        out.extend(round_.product_segments())
    return out


def _base_segments(kind: str) -> list[Segment]:
    parts = [Segment("shared ligation anchor", SHARED_ANCHOR, "r2")]
    if kind == "RNA":
        parts.extend([
            Segment("mosaic-end-derived RT handle", ME, "me"),
            Segment("UMI", "U" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
            Segment("internal biotin-dT", "T"),
            Segment("oligo-dT and VN anchor", RT_DT, placeholder=True),
            Segment("first-strand cDNA", "X" * 28, placeholder=True),
        ])
    else:
        parts.extend([
            Segment("mosaic end", ME, "me"),
            Segment("accessible genomic DNA", "X" * 36, placeholder=True),
        ])
    return parts


def _molecule_scene(kind: str, rounds: int) -> Scene:
    parts = [*_prefix(rounds), *_base_segments(kind)]
    sc = Scene()
    strand = kind.lower()
    sc.strand(strand, parts, label=f"{kind} ligatable strand",
              mod5="/5Phos/" if rounds < 3 else "")
    if kind == "RNA":
        sc.mark(strand, "internal biotin-dT", "internal biotin")
    previous = "shared ligation anchor"
    for round_ in ROUNDS[:rounds]:
        sc.junction(strand, f"{round_.name} ligating arm", previous,
                    f"{round_.name} ligation")
        previous = f"{round_.name} outgoing arm"
    sc.labels(strand)
    return sc


def _hybridized_scene(kind: str, rounds: int) -> Scene:
    """Pre-ligation barcode chain: every backbone oligo remains a separate strand."""
    sc = Scene()
    substrate = f"{kind.lower()} substrate"
    sc.strand(substrate, _base_segments(kind), label="",
              rev=True, mod5="/5Phos/")
    prior_strand = substrate
    prior_arm = "shared ligation anchor"
    for round_ in ROUNDS[:rounds]:
        linker_name = f"{kind.lower()} {round_.name} linker"
        barcode_name = f"{kind.lower()} {round_.name} barcode"
        linker = [
            Segment(f"{round_.name} prior-molecule splint",
                    round_.linker[:len(round_.previous_5_arm)]),
            Segment(f"{round_.name} barcode splint",
                    round_.linker[len(round_.previous_5_arm):]),
        ]
        sc.anneal(linker_name, linker, to=prior_strand,
                  pair=(f"{round_.name} prior-molecule splint", prior_arm),
                  label="")
        sc.anneal(barcode_name, round_.product_segments(), to=linker_name,
                  pair=(f"{round_.name} ligating arm",
                        f"{round_.name} barcode splint"),
                  label="",
                  mod5="/5Phos/" if round_.phosphorylated_5 else "")
        prior_strand = barcode_name
        prior_arm = f"{round_.name} outgoing arm"
    sc.note(substrate, "annealed barcode chain; backbone nicks remain unsealed")
    sc.labels(substrate)
    for round_ in ROUNDS[:rounds]:
        sc.labels(f"{kind.lower()} {round_.name} barcode")
    return sc


def _both_modalities(rounds: int):
    return [*_molecule_scene("RNA", rounds).rows(),
            *_molecule_scene("ATAC", rounds).rows()]


def workflow():
    wf = Workflow(MolecularState("Fixed permeabilized cell", tuple(_initial_rows())))
    wf.react(
        "Tagment accessible chromatin with the ligatable Nextera end",
        _molecule_scene("ATAC", 0).rows(), name="Tn5-tagged chromatin",
        note="The selected Tn5 product exposes the same 15-nt ligation anchor used by the RT primer; other Tn5 end combinations do not both ligate and amplify.")
    wf.react(
        "Reverse-transcribe poly(A) RNA with the phosphorylated biotin-UMI primer",
        _both_modalities(0), name="RNA and ATAC molecules with a shared ligation anchor",
        note="The exact RT primer carries a 10-nt UMI, internal biotin-dT and oligo-dT/VN tract. Its 5′ constant arm matches the ligatable Tn5 end.")
    for number, round_ in enumerate(ROUNDS, 1):
        wf.react(
            f"Hybridize {round_.name}, block free adapters, pool and resplit",
            [*_hybridized_scene("RNA", number).rows(),
             *_hybridized_scene("ATAC", number).rows()],
            name=f"After {round_.name} hybridization",
            note="The linker is computed from the reverse complements of the incoming molecule arm and barcode 3′ arm. A complementary blocking oligo sequesters adapters that did not bind a molecule before cells are pooled.")
    wf.react(
        "Add T4 DNA ligase once to seal all three barcode junctions",
        _both_modalities(3), name="Covalently barcoded RNA and ATAC molecules",
        note="The three marked bonds are created only after the serial hybridization and blocking rounds have assembled the complete chain.")
    wf.react(
        "Release molecules, reverse crosslinks and separate the two modalities",
        modality_split("streptavidin-captured RNA-derived cDNA",
                       "ATAC-derived genomic DNA"),
        name="Separated RNA and ATAC fractions",
        note="Internal biotin on the RT primer enables physical separation before modality-specific amplification. Template switching and RNA-library PCR follow on the captured branch.")
    return wf
