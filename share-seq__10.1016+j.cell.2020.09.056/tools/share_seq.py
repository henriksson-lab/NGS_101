"""SHARE-seq joint RNA/accessibility combinatorial indexing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import MolecularState, Row, Scene, Segment, Workflow, feature
from multimodal_spatial import (HybridizationLigationRound, finish_pcr_branch,
                                modality_split)
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE="SHARE-seq"
NOTES="01_share-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma et al., <i>Cell</i> (2020)</a>.'
SUMMARY="Fixed permeabilized cells are first tagmented and reverse-transcribed. Three split-pool hybridization rounds assemble the same barcode chain on ATAC fragments and biotinylated cDNA; blocking oligos sequester free adapters, then one ligation seals the chain before the modalities are separated."
CAVEAT=("The constant arms, linkers and ordering modifications come from the published "
        "Supplementary Table S1. Its 96 well-specific 8-nt values per round are represented "
        "as placeholders in the diagram. <br><b>Source-version discrepancy:</b> the "
        "<a href='https://doi.org/10.1101/2020.06.17.156943'>June 2020 preprint</a> "
        "Methods, section (1), calls for the <i>round 2</i> linker in the round 1 plate "
        "(PDF lines 722–725) and describes the incoming molecular end as a 3′ overhang "
        "(lines 718–721). The <a href='https://doi.org/10.1016/j.cell.2020.09.056'>"
        "published STAR Methods</a>, ‘Preparing oligonucleotides for ligations’, calls for "
        "the <i>round 1</i> linker and a 5′ overhang (JATS line 213); Table S1 ‘bulk oligos’ "
        "lines 3–5 prints distinct linkers for all three rounds. The model uses the "
        "published linker and oligo-table geometry. The preprint's intended wording "
        "cannot be independently established from its methods text alone.")
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


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","fixed permeabilized cells containing chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=_initial_rows(),primary=True)
    transposome=g.input("transposome","ligatable-anchor-loaded Nextera transposome",MaterialKind.ADAPTER)
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","cells retaining Tn5-tagged chromatin and RNA",MaterialKind.MOLECULAR,state=[*_molecule_scene("ATAC",0).rows(),Row(chunks=[("poly(A) RNA remains unindexed in the intact cell",None,False)])])
    twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tagment accessible chromatin with the ligatable Nextera end",StepKind.TAGMENTATION,((cells,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not expose the selected ligation anchor on an amplifiable fragment")
    rtprimer=g.input("rt_primer","phosphorylated biotin-UMI oligo(dT)/VN RT primer",MaterialKind.OLIGO)
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    indexed=g.product("indexed","RNA and ATAC molecules exposing the shared ligation anchor",MaterialKind.MOLECULAR,state=_both_modalities(0))
    g.react("reverse_transcribe","Reverse-transcribe poly(A) RNA with the shared-anchor primer",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(indexed,))
    current=indexed
    reviewed=["cells","transposome","tn5","tagment","rt_primer","rt","reverse_transcribe"]
    for number,round_ in enumerate(ROUNDS,1):
        barcode=g.input(f"barcode_{number}",f"{round_.name} barcode oligo",MaterialKind.OLIGO)
        linker=g.input(f"linker_{number}",f"{round_.name} splint linker",MaterialKind.OLIGO)
        hybrid=g.product(f"hybrid_{number}",f"RNA and ATAC with {round_.name} annealed",MaterialKind.MOLECULAR,state=[*_hybridized_scene("RNA",number).rows(),*_hybridized_scene("ATAC",number).rows()])
        g.react(f"hybridize_{number}",f"Hybridize {round_.name} barcode and splint",StepKind.HYBRIDIZATION,((current,InputRole.TEMPLATE),(barcode,InputRole.STRUCTURAL),(linker,InputRole.STRUCTURAL)),(hybrid,))
        blocker=g.input(f"blocker_{number}",f"{round_.name} free-adapter blocking oligo",MaterialKind.OLIGO)
        cleaned=g.product(f"cleaned_{number}",f"{round_.name} assemblies after free-adapter blocking",MaterialKind.MOLECULAR,state=hybrid.state)
        blocked=g.product(f"blocked_{number}",f"blocked unbound {round_.name} assemblies",MaterialKind.WASTE)
        g.react(f"block_{number}",f"Sequester unbound {round_.name} adapters",StepKind.CLEANUP,((hybrid,InputRole.TEMPLATE),(blocker,InputRole.STRUCTURAL)),(cleaned,blocked),changes_molecule=False)
        g.discard(blocked,"is an unbound adapter assembly sequestered before repooling")
        redistribution=g.input(f"redistribution_{number}",f"{round_.name} pool-and-resplit condition",MaterialKind.CONDITION)
        redistributed=g.product(f"redistributed_{number}",f"molecules redistributed after {round_.name}",MaterialKind.MOLECULAR,state=cleaned.state)
        excluded=g.product(f"excluded_{number}",f"unresolved {round_.name} redistribution material",MaterialKind.WASTE)
        g.react(f"redistribute_{number}",f"Pool and resplit after {round_.name}",StepKind.PARTITION,((cleaned,InputRole.TEMPLATE),(redistribution,InputRole.CONDITION)),(redistributed,excluded),changes_molecule=False)
        g.discard(excluded,"does not retain an interpretable ordered well path")
        reviewed.extend((barcode.id,linker.id,f"hybridize_{number}",blocker.id,f"block_{number}",redistribution.id,f"redistribute_{number}"))
        current=redistributed
    ligase=g.input("ligase","T4 DNA ligase",MaterialKind.CATALYST)
    sealed=g.product("sealed","RNA and ATAC molecules with three sealed barcode junctions",MaterialKind.MOLECULAR,state=_both_modalities(3))
    g.react("seal","Seal all three preassembled barcode-chain nicks",StepKind.NICK_SEALING,((current,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(sealed,))
    release=g.input("release_condition","lysis and crosslink-reversal condition",MaterialKind.CONDITION)
    released=g.product("released","released covalently barcoded RNA and ATAC molecules",MaterialKind.MOLECULAR,state=sealed.state)
    g.react("release","Lyse cells and reverse crosslinks",StepKind.EXTRACTION,((sealed,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(released,))
    streptavidin=g.input("streptavidin","streptavidin selection matrix",MaterialKind.SELECTION)
    rna=g.product("rna_branch","biotin-captured RNA-derived cDNA",MaterialKind.MOLECULAR,state=_molecule_scene("RNA",3))
    atac=g.product("atac_branch","ATAC-derived genomic-DNA flow-through",MaterialKind.MOLECULAR,state=_molecule_scene("ATAC",3))
    g.react("separate","Separate RNA-derived and ATAC-derived molecules",StepKind.AFFINITY_CAPTURE,((released,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(rna,atac),changes_molecule=False)
    tso=g.input("tso","SHARE-seq template-switching oligo",MaterialKind.OLIGO)
    switched=g.product("switched","template-switched RNA-derived cDNA",MaterialKind.MOLECULAR,state=rna.state)
    g.react("template_switch","Template-switch the captured RNA branch",StepKind.TEMPLATE_SWITCHING,((rna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    reviewed.extend(("ligase","seal","release_condition","release","streptavidin","separate","tso","template_switch"))
    finish_pcr_branch(g,"rna",switched,RNA,RNA_P,reviewed,label="SHARE-seq RNA")
    finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="SHARE-seq ATAC")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
