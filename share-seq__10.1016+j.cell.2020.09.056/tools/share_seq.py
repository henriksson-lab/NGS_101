"""SHARE-seq joint RNA/accessibility combinatorial indexing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from dataclasses import dataclass
from protocol_model import role
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal

TITLE="SHARE-seq"
NOTES="01_share-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma et al., <i>Cell</i> (2020)</a>.'
SUMMARY="Fixed permeabilized cells are first tagmented and reverse-transcribed. Three split-pool hybridization rounds assemble the same barcode chain on ATAC fragments and biotinylated cDNA; blocking oligos sequester free adapters, then one ligation seals the chain before the modalities are separated."
CAVEAT=("Published STAR Methods and TableS1 establish the shared5′ligation anchor, distinct three-round linkers and final Ad1/P7 ends. Superseded preprint wording is not a published-source conflict. Native first-cycle/TSO/end-gap proofs remain unassigned; source-primer-derived final libraries are separate references. Blocking is hybridization, not physical disposal; well-path/collision filtering remains computational. Source NextSeq30/99/8/30 is selected; NovaSeq50/99/8/50 is separate. Sequencing primer compatibility is inferred, while acquired cycle allocation is source-positive.")
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

@dataclass(frozen=True)
class SourceRound:
    name: str
    previous_5_arm: str
    barcode_5_arm: str
    barcode_3_arm: str
    linker: str
    feature: object
    phosphorylated_5: bool
    barcode_length: int
    blocking_oligo: str

    def __post_init__(self):
        if self.linker != revcomp(self.previous_5_arm)+revcomp(self.barcode_3_arm):
            raise ValueError("TableS1 linker does not bridge ordered source arms")
        if self.blocking_oligo != revcomp(self.linker):
            raise ValueError("printed blocker does not complement source linker")

    def product_segments(self):
        return [Segment(self.name+" outgoing arm",self.barcode_5_arm),Segment(self.name+" barcode","N"*self.barcode_length,"cbc",placeholder=True,feature=self.feature),Segment(self.name+" ligating arm",self.barcode_3_arm)]

ROUNDS = (
    SourceRound(
        "round 1", SHARED_ANCHOR, "CGCGCTGCATACTTG", "CCCATGATCGTCCGA",
        "CCGAGCCCACGAGACTCGGACGATCATGGG", BC1, phosphorylated_5=True, barcode_length=8, blocking_oligo="CCCATGATCGTCCGAGTCTCGTGGGCTCGG"),
    SourceRound(
        "round 2", "CGCGCTGCATACTTG", "CATCGGCGTACGACT", "ATCCACGTGCTTGAG",
        "CAAGTATGCAGCGCGCTCAAGCACGTGGAT", BC2, phosphorylated_5=True, barcode_length=8, blocking_oligo="ATCCACGTGCTTGAGCGCGCTGCATACTTG"),
    SourceRound(
        "round 3", "CATCGGCGTACGACT", "CAAGCAGAAGACGGCATACGAGAT",
        "GTGGCCGATGTTTCG", "AGTCGTACGCCGATGCGAAACATCGGCCAC", BC3, phosphorylated_5=False, barcode_length=8, blocking_oligo="GTGGCCGATGTTTCGCATCGGCGTACGACT"),
)

# Published STAR Methods "scATAC-seq library preparation", "Tagmentation and scRNA-seq
# library preparation" and "Quantification and sequencing" (JATS lines 241-249) with
# Table S1: both libraries are amplified with a sample-barcoded Ad1 primer
# (P5-[8-nt i5]-TCGTCGGCAGCGTCAGATGTGTAT, sheet "Ad1") and the P7 primer, which is the
# round 3 barcode oligo's 5' arm.  The three hybridization barcodes therefore lie
# between the shared S7-derived anchor and P7, and the 99-cycle Index 1 read covers
# them; there is no separate i7 index.  (Model correction 2026-10-17: earlier versions
# drew generic TruSeq/Nextera dual-index libraries with the barcodes inside the insert.)
AD1_I5_FEATURE = feature("ad1_sublibrary_index", "sample_index", "whitelist",
                         whitelist="Table S1 sheet Ad1, v2_Ad1.01-96",
                         note="PCR sub-library barcode; the paper combines it with the three hybridization barcodes")


def _p7_barcode_chain() -> list[Segment]:
    """Barcode chain read 5'->3' on the P5-oriented strand: reverse complement of the ligated strand."""
    out: list[Segment] = [Segment("shared anchor reverse complement (S7 part)",
                                  revcomp(SHARED_ANCHOR), "r2")]
    for round_ in ROUNDS:
        out.append(Segment(f"{round_.name} ligating arm reverse complement",
                           revcomp(round_.barcode_3_arm)))
        out.append(Segment(f"{round_.name} barcode", "B" * round_.barcode_length, "cbc",
                           placeholder=True, feature=round_.feature))
        out.append(Segment(f"{round_.name} outgoing arm reverse complement",
                           revcomp(round_.barcode_5_arm),
                           "p7" if round_ is ROUNDS[-1] else None))
    return out


def _final_library(name: str, insert: list[Segment]) -> Construct:
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Ad1 sub-library barcode (i5)", "J" * 8, "cbc", placeholder=True,
                feature=AD1_I5_FEATURE),
        Segment("Read 1 arm (S5)", nx.S5, "s5"),
        Segment("mosaic end", nx.ME, "me"),
        *insert,
        Segment("mosaic end reverse complement", nx.ME_RC, "me"),
        *_p7_barcode_chain(),
    ], name=name)


# Library PCR primers as printed in Table S1 (sheet "Ad1" family and "bulk oligos" P7).
AD1_PRIMER = Construct([
    Segment("P5", il.P5, "p5"),
    Segment("Ad1 sub-library barcode (i5)", "J" * 8, "cbc", placeholder=True,
            feature=AD1_I5_FEATURE),
    Segment("Read 1 arm (S5) and mosaic-end start", nx.S5 + nx.ME[:10], "s5"),
], name="v2_Ad1 primer family")
P7_PRIMER = Construct([Segment("P7", il.P7, "p7")], name="P7 primer")
LIBRARY_PCR_SOURCE = ("Ma et al. (2020), STAR Methods (JATS lines 241 and 249) and "
                      "Table S1 (sheets 'Ad1' and 'bulk oligos')")


RNA = _final_library("SHARE-seq RNA library", [
    Segment("Tn5-fragmented cDNA (mRNA sense)", "X" * 30, placeholder=True),
    Segment("VN complement", "NB", placeholder=True),
    Segment("poly(A) complement of dT tract and biotin-dT", "A" * 15),
    Segment("UMI", "U" * 10, "umi", placeholder=True, feature=UMI_FEATURE),
])
ATAC = _final_library("SHARE-seq ATAC library", [
    Segment("accessible DNA", "X" * 34, placeholder=True),
])
# Primers are not named in the paper; the standard Nextera set matches the Ad1/S7-ME
# adapter sequences.  NextSeq 150-cycle kit configuration (the NovaSeq S1 alternative
# is 50/99/8/50).
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
READ_LENGTHS_NEXTSEQ = {"Read 1": 30, "Index 1 (i7)": 99, "Index 2 (i5)": 8, "Read 2": 30}
RNA_P = ATAC_P = SEQ_PRIMERS
_READ_NOTE = (" Published run: NextSeq 150-cycle High-Output (Read 1 30, Index 1 99, Index 2 8,"
              " Read 2 30) or NovaSeq S1 200-cycle (50/99/8/50); the NextSeq lengths are"
              " drawn. Index 1 reads the three 8-nt hybridization barcodes with their"
              " 15-nt arms (99 cycles); Index 2 reads the Ad1 sub-library barcode. The paper"
              " does not name the sequencing primers; the standard Nextera primers matching"
              " these adapters are assumed.")
FINAL_LIBRARIES=(
    ("Final RNA library",RNA,RNA_P,"Ad1/P7 PCR after Read 1-only Tn5 tagmentation of amplified cDNA. Read 2 starts with the 10-nt UMI followed by the dT tract (the pipeline requires TTTTTT at Read 2 bases 11-16)."+_READ_NOTE,"The RNA library's declared primer sites are computed on the completed duplex."),
    ("Final ATAC library",ATAC,ATAC_P,"Ad1/P7 PCR of the streptavidin supernatant; the same three hybridization barcodes identify the cell of origin of each accessible-DNA fragment."+_READ_NOTE,"The ATAC library's declared primer sites are computed on the completed duplex."),)

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
        f"<b>P7 primer</b> — 5′-{il.P7}-3′ (Table S1 bulk oligos; identical to the round 3 barcode 5′ arm)",
        f"<b>Ad1 primer family</b> — 5′-{il.P5}[8-nt sub-library barcode]{nx.S5 + nx.ME[:10]}-3′ (Table S1 sheet Ad1, 96 values)",
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


def _native_rna_rows():
    sc=Scene();sc.strand("nativeRNA",[role("independent native RNA template; recovery/degradation unassigned")],label="RNA is separate from its DNA copy");return sc.rows()


def _captured_rna():
    sc=_molecule_scene("RNA",3);sc.strand("residual RNA",[role("independent RNA template; exact retention through affinity capture/washes unassigned")],label="possible carried RNA, not a covalent cDNA segment");return sc


def _both_modalities(rounds: int):
    return [*_molecule_scene("RNA", rounds).rows(),
            *_molecule_scene("ATAC", rounds).rows(),*_native_rna_rows()]


TSO_DNA = RNA_PCR + "GAAT"   # DNA part of the TSO; the 3' end is rG rG +G (written above)


def _amplified_cdna() -> Construct:
    """Ligated first strand plus TSO complement, drawn in the ligated-strand orientation."""
    return Construct([
        *_prefix(3), *_base_segments("RNA"),
        Segment("complement of TSO (CCC opposite rGrG+G)", revcomp(TSO_DNA + "GGG")),
    ], name="amplified SHARE-seq cDNA")


def _tagmented_cdna() -> Scene:
    """Native partial-strand input; never truncate the expected completed library."""
    sc=Scene()
    sc.strand("barcode bearing strand",[*_prefix(3),*_base_segments("RNA")],label="barcode-bearing cDNA fragment; native break extent unassigned")
    sc.strand("Read1 bearing strand",[Segment("transferred Read1 oligo",READ1),role("opposite native cDNA strand; paired span/gaps unassigned")],label="separate donor-bearing strand; no manufactured completed duplex")
    return sc


def _loaded(read1_only=False):
 sc=Scene()
 for key,seq,phos in (("Read1",READ1,False),) if read1_only else (("Read1",READ1,False),("Read2",RT_HANDLE,True)):
  # Published blocked14nt ME complement pairs to the final14 bases, with5′P, PT bonds and3′ddC.
  sc.strand(key,[Segment(key+" prefix",seq[:-14]),Segment(key+" ME paired end",seq[-14:])],label="source transferred oligo",mod5="p" if phos else "")
  sc.anneal(key+" complement",[Segment(key+" blocked ME complement","CTGTCTCTTATACA",note="Source phosphorothioates between the printed starred residues")],key,pair=(key+" blocked ME complement",key+" ME paired end"),label="source loading strand",mod5="p",mod3="ddC")
 return sc


def _rt_primer():
 sc=Scene();sc.strand("RT",_base_segments("RNA")[:-1],label="source RT oligo",mod5="p");sc.mark("RT","internal biotin-dT","internal biotin-dT");return sc


def _ss_switched():
 sc=_molecule_scene("RNA",3)
 sc.strands["rna"].segs.append(Segment("copied TSO boundary",revcomp(TSO_DNA+"GGG"),inferred=True))
 sc.strand("RNA template",[role("native RNA template retained; removal/degradation unassigned")],label="independent RNA, not the copied strand")
 return sc


def _clean(g,key,template,label,state=None):
 selection=g.input(key+"_selection",label,MaterialKind.SELECTION)
 kept=g.product(key+"_retained","recovered "+template.label,MaterialKind.MOLECULAR,state=template.state if state is None else state)
 waste=g.product(key+"_waste","physical supernatant/washes; composition unassigned",MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,"Physical non-retained fraction; class-specific removal is not inferred.")
 return kept


def _controlled_pcr(g,key,template,f,r,enzyme,state,qc_volume,qc_dye):
 residual=g.product(key+"_residual","co-carried RNA/other polymers and off-target products; physical recovery unassigned",MaterialKind.MOLECULAR,state=template.state)
 amp=g.product(key+"_five","5-cycle source PCR material",MaterialKind.MOLECULAR,state=state)
 g.react(key+"_five_cycles","Initial5cycles source "+key+" PCR",StepKind.PCR,((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(amp,residual),geometry=PCRGeometry.unavailable(template.state,reason="Actual native precursor; complete first-cycle/end geometry remains unassigned."))
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"PCR specificity does not physically remove co-carried RNA or off-target polymers; recovery through later cleanups unassigned."))
 production=g.product(key+"_production","remaining5-cycle material for production PCR",MaterialKind.MOLECULAR,state=state)
 qc=g.product(key+"_qc_aliquot",qc_volume+" analytical qPCR aliquot",MaterialKind.MOLECULAR,state=state)
 g.react(key+"_qc_split","Take source qPCR aliquot after5cycles",StepKind.SPLIT,((amp,InputRole.TEMPLATE),),(production,qc),changes_molecule=False)
 dye=g.input(key+"_qc_dye",qc_dye+" qPCR cocktail;20cycles to saturation",MaterialKind.CONDITION)
 qc_product=g.product(key+"_qc_product","analytical qPCR product measuring1/3saturation cycle requirement",MaterialKind.QC,state=state)
 g.react(key+"_qpcr","Measure remaining cycle requirement by qPCR",StepKind.PCR,((qc,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST),(dye,InputRole.CONDITION)),(qc_product,),geometry=PCRGeometry.unavailable(qc.state,reason="Actual analytical aliquot; qPCR product geometry not separately proven."))
 g.mark_terminal(Terminal(qc_product.id,Disposition.QC,"Analytical amplification used to choose1/3saturation cycle count, not the production library."))
 extra=g.input(key+"_extra_cycles","extra cycle count selected from qPCR1/3saturation; sample-specific",MaterialKind.CONDITION)
 completed=g.product(key+"_amplified","source qPCR-controlled production amplicons",MaterialKind.MOLECULAR,state=state)
 g.react(key+"_extra_pcr","Complete production PCR with qPCR-defined extra cycles",StepKind.PCR,((production,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST),(extra,InputRole.CONDITION)),(completed,),geometry=PCRGeometry.unavailable(production.state,reason="Actual5-cycle production aliquot; native selected ends remain bounded."))
 return completed


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","fixed permeabilized cells/nuclei, selected cell-line workflow",MaterialKind.BIOLOGICAL,state=_initial_rows(),primary=True)
 loaded=g.input("transposome","source Read1 and5′P Read2 with5′P/PT/3′ddC blocked ME complements",MaterialKind.ADAPTER,state=_loaded())
 tn5=g.input("tn5","in-houseTn5,37°C30min",MaterialKind.CATALYST)
 tagged=g.product("tagged","Tn5-tagged chromatin with native RNA retained",MaterialKind.MOLECULAR,state=[*_molecule_scene("ATAC",0).rows(),*_native_rna_rows()])
 other_ends=g.product("other_ends","other Tn5 end classes",MaterialKind.MOLECULAR)
 g.react("tagment","Tagment fixed chromatin with source transferred strands",StepKind.TAGMENTATION,((cells,InputRole.TEMPLATE),(loaded,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other_ends))
 g.mark_terminal(Terminal(other_ends.id,Disposition.UNRESOLVED,"End eligibility does not establish zero amplification or physical exclusion; recovery in subsequent nuclear washes remains unassigned."))
 tagged=_clean(g,"wash_tagment",tagged,"NIB pellet/wash after transposition")
 primer=g.input("rt_primer","source5′P anchor/N10UMI/biotin-dT/T14VN",MaterialKind.OLIGO,state=_rt_primer())
 rt=g.input("rt","MaximaHMinus RT;50°C then3temperature-ramp cycles then50°C5min",MaterialKind.CATALYST)
 indexed=g.product("indexed","ss-cDNA and independent ATAC populations with shared anchor",MaterialKind.MOLECULAR,state=_both_modalities(0))
 g.react("reverse_transcribe","RT native poly(A)RNA in fixed cells",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(indexed,))
 current=_clean(g,"wash_rt",indexed,"NIB nuclear pellet/wash after RT")
 for n,round_ in enumerate(ROUNDS,1):
  sc=Scene();sc.strand("barcode",round_.product_segments(),label=round_.name,mod5="p" if round_.phosphorylated_5 else "")
  barcode=g.input("barcode_"+str(n),round_.name+" source N8 barcode",MaterialKind.OLIGO,state=sc)
  linker=g.input("linker_"+str(n),round_.name+" source splint",MaterialKind.OLIGO,state=round_.linker)
  hybrid=g.product("hybrid_"+str(n),"independent RNA/ATAC pre-ligation assemblies",MaterialKind.MOLECULAR,state=[*_hybridized_scene("RNA",n).rows(),*_hybridized_scene("ATAC",n).rows(),*_native_rna_rows()])
  g.react("hybridize_"+str(n),"Split among96wells and hybridize "+round_.name,StepKind.HYBRIDIZATION,((current,InputRole.TEMPLATE),(barcode,InputRole.STRUCTURAL),(linker,InputRole.STRUCTURAL)),(hybrid,))
  blocker=g.input("blocker_"+str(n),round_.name+" source complementary blocking oligo",MaterialKind.OLIGO,state=round_.blocking_oligo)
  blocked=g.product("blocked_"+str(n),"cell-associated assemblies plus blocked free oligos; nicks unsealed",MaterialKind.MOLECULAR,state=hybrid.state)
  g.react("block_"+str(n),"Hybridize blocking strand then pool all96wells",StepKind.HYBRIDIZATION,((hybrid,InputRole.TEMPLATE),(blocker,InputRole.STRUCTURAL)),(blocked,))
  current=blocked
 current=_clean(g,"wash_hybridizations",current,"After round3 pool, pellet and wash twice with NIB")
 ligase=g.input("ligase","20U/µlT4DNA ligase;RT30min",MaterialKind.CATALYST)
 sealed=g.product("sealed","covalently sealed three-round cDNA/ATAC barcode chains",MaterialKind.MOLECULAR,state=_both_modalities(3))
 g.react("seal","Seal3preassembled barcode-chain nicks once",StepKind.NICK_SEALING,((current,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(sealed,))
 sealed=_clean(g,"wash_ligation",sealed,"NIB wash after ligation")
 aliquot=g.input("cell_aliquot","1000–20000barcoded cells per tube; other aliquots outside selected route",MaterialKind.CONDITION)
 selected=g.product("selected_cells","selected barcoded cell aliquot",MaterialKind.MOLECULAR,state=sealed.state)
 remainder=g.product("other_cells","other barcoded aliquots",MaterialKind.MOLECULAR,state=sealed.state)
 g.react("aliquot_cells","Aliquot sealed barcoded cells",StepKind.SAMPLING,((sealed,InputRole.TEMPLATE),(aliquot,InputRole.CONDITION)),(selected,remainder),changes_molecule=False)
 g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"Other library aliquots are outside this representation; no computed path exclusion is physical discard."))
 release=g.input("release_condition","proteinaseK/SDS55°C1h;PMSFRT10min",MaterialKind.CONDITION)
 released=g.product("released","released barcode-bearing cDNA/ATAC populations",MaterialKind.MOLECULAR,state=selected.state)
 g.react("release","Reverse crosslinks and inactivate proteinaseK",StepKind.EXTRACTION,((selected,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(released,))
 beads=g.input("streptavidin","MyOneC1 streptavidin beads;RT60min",MaterialKind.SELECTION)
 rna=g.product("rna_branch","biotin-captured ss-cDNA/RNA reaction material",MaterialKind.MOLECULAR,state=_captured_rna())
 atac=g.product("atac_branch","supernatant containing genomic ATAC fragments",MaterialKind.MOLECULAR,state=_molecule_scene("ATAC",3))
 g.react("affinity_split","Physically retain biotinylated cDNA on beads; ATAC in supernatant",StepKind.AFFINITY_CAPTURE,((released,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,atac),changes_molecule=False)
 atac=_clean(g,"atac_minelute",atac,"QIAGENMinElute PCRcleanup;20µlTris elution")
 rna=_clean(g,"rna_bead_wash",rna,"ThreeB&W-T bead washes plus oneSTE wash")
 tso_sc=Scene();tso_sc.strand("TSO",[Segment("TSO DNA",TSO_DNA),Segment("terminalrGrG+G","GGG",note="first2G RNA; lastG LNA")],label="source TSO modified tail")
 tso=g.input("tso","TableS1 TSO withrGrG+G",MaterialKind.OLIGO,state=tso_sc)
 switched=g.product("switched","actual switched ss-cDNA beside retained RNA",MaterialKind.MOLECULAR,state=_ss_switched(),inferred=True)
 g.react("template_switch","INFERRED copied-end reference — on-beadMaxima template switch:RT30min,42°C90min",StepKind.TEMPLATE_SWITCHING,((rna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 switched=_clean(g,"wash_switch",switched,"Magnetic supernatant removal andSTE wash after template switching")
 p7=g.input("p7","sourceP7 primer",MaterialKind.OLIGO,state=P7_PRIMER)
 rnap=g.input("rna_pcr_primer","sourceRNA_PCR_primer",MaterialKind.OLIGO,state=RNA_PCR)
 kapa=g.input("kapa","KAPAHifi cDNA PCRmix",MaterialKind.CATALYST)
 amp=_controlled_pcr(g,"cdna",switched,p7,rnap,kapa,_amplified_cdna(),"2.5µl","1×EvaGreen")
 amp=_clean(g,"cdna_ampure",amp,"Selected cell-line0.8×AMPure;10µlTris elution (primary-tissue0.6×separate)")
 amount=g.input("cdna_amount","50ng amplified cDNA for fragmentation",MaterialKind.CONDITION)
 input_cdna=g.product("cdna_input","50ng cDNA aliquot",MaterialKind.MOLECULAR,state=amp.state)
 rest=g.product("cdna_remainder","other recovered cDNA",MaterialKind.MOLECULAR,state=amp.state)
 g.react("sample_cdna","Take50ng cDNA",StepKind.SAMPLING,((amp,InputRole.TEMPLATE),(amount,InputRole.CONDITION)),(input_cdna,rest),changes_molecule=False)
 g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,"cDNA outside selected50ng input; fate not assigned."))
 read1=g.input("rna_transposome","sourceRead1-onlyTn5 with blockedME complement",MaterialKind.ADAPTER,state=_loaded(True))
 rna_tn5=g.input("rna_tn5","Read1-onlyTn5;55°C5min",MaterialKind.CATALYST)
 fragments=g.product("rna_fragments","native Read1-tagged cDNA reference with barcode end; gap extent unassigned",MaterialKind.MOLECULAR,state=_tagmented_cdna())
 off=g.product("rna_off_targets","other cDNA fragments/end classes",MaterialKind.MOLECULAR)
 g.react("rna_tagment","Tagment cDNA withRead1-onlytransposome",StepKind.TAGMENTATION,((input_cdna,InputRole.TEMPLATE),(read1,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(fragments,off))
 g.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,"End eligibility is not physical removal or proof of zero amplification; recovery through Zymo/AMPure unassigned."))
 fragments=_clean(g,"rna_zymo",fragments,"ZymoDNAclean&Concentrator;10µlTris elution")
 ad1=g.input("ad1","sourceAd1sample-indexed primer familyN8",MaterialKind.OLIGO,state=AD1_PRIMER)
 pol=g.input("library_polymerase","NEBNextHighFidelity2×PCRMasterMix",MaterialKind.CATALYST)
 amplified=g.product("rna_amplified","source RNA indexed amplicons",MaterialKind.MOLECULAR,state=RNA)
 g.react("rna_index_pcr","72°C5min then7cyclesAd1/P7PCR",StepKind.PCR,((fragments,InputRole.TEMPLATE),(ad1,InputRole.STRUCTURAL),(p7,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(fragments.state,reason="Actual nativeRead1-tagged input; precise gap/first-cycle derivation unassigned."))
 rna_clean=_clean(g,"rna_final_cleanup",amplified,"0.7×AMPure;10µlTris elution")
 rna_library=g.product("rna_library",RNA.name,MaterialKind.LIBRARY,state=RNA)
 g.react("rna_qc","Quantify recovered RNA library withKAPA;retain for sequencing",StepKind.GENERIC,((rna_clean,InputRole.TEMPLATE),),(rna_library,),changes_molecule=False)
 atac_amp=_controlled_pcr(g,"atac",atac,ad1,p7,pol,ATAC,"25µl","0.6×SYBRgreen")
 atac_library=g.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC)
 g.react("atac_qc","Quantify source ATAC library withKAPA;retain for sequencing",StepKind.GENERIC,((atac_amp,InputRole.TEMPLATE),),(atac_library,),changes_molecule=False)
 for lib,endpoint,primers in ((rna_library,RNA,RNA_P),(atac_library,ATAC,ATAC_P)):
  g.sequence(lib,SequencingPlan(library=endpoint,primers=tuple(primers),acquired_roles=tuple(READ_LENGTHS_NEXTSEQ),read_lengths=READ_LENGTHS_NEXTSEQ,limitation="Compatible Nextera sequencing primers inferred from source arms; acquired NextSeq30/99/8/30 source-positive, NovaSeq50/99/8/50 separate."))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=("cells","transposome","tn5","tagment","wash_tagment_selection","wash_tagment","rt_primer","rt","reverse_transcribe","wash_rt_selection","wash_rt","barcode_1","linker_1","hybridize_1","blocker_1","block_1","barcode_2","linker_2","hybridize_2","blocker_2","block_2","barcode_3","linker_3","hybridize_3","blocker_3","block_3","wash_hybridizations_selection","wash_hybridizations","ligase","seal","wash_ligation_selection","wash_ligation","cell_aliquot","aliquot_cells","release_condition","release","streptavidin","affinity_split","atac_minelute_selection","atac_minelute","rna_bead_wash_selection","rna_bead_wash","tso","template_switch","wash_switch_selection","wash_switch","p7","rna_pcr_primer","kapa","cdna_five_cycles","cdna_qc_split","cdna_qc_dye","cdna_qpcr","cdna_extra_cycles","cdna_extra_pcr","cdna_ampure_selection","cdna_ampure","cdna_amount","sample_cdna","rna_transposome","rna_tn5","rna_tagment","rna_zymo_selection","rna_zymo","ad1","library_polymerase","rna_index_pcr","rna_final_cleanup_selection","rna_final_cleanup","rna_qc","rna_library","atac_five_cycles","atac_qc_split","atac_qc_dye","atac_qpcr","atac_extra_cycles","atac_extra_pcr","atac_qc","atac_library")
