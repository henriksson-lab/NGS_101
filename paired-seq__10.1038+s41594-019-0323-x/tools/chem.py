"""Paired-seq joint RNA and chromatin accessibility.

Exact oligos are from Zhu et al. Supplementary Tables 1 and 2. The molecular
model uses the zero-phase members of the round-4 plate; other wells insert one
to three phase bases between BC4 and its ligating arm.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      duplex_rows, feature, revcomp)
from dataclasses import dataclass,replace
from protocol_model import molecule,role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind, SequencingPlan, Disposition, Terminal, Highlight

SOURCE_OLIGO_EXCEPTIONS={"bc1_transposome":"Published SupplementaryTable1 writes pMENTs terminalddC as a nucleotide; scraper._letters strips it as a modification without its physicalC. Local source test verifies whole modified pMENTs plus source barcoded transferred strand directly against original Tables1/2. Exact physical terminalC remains; no alternative sequence supplied."}

TITLE = "Paired-seq — combinatorial RNA and chromatin accessibility"
NOTES = "01_paired-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41594-019-0323-x">'
          'Zhu et al., <i>Nature Structural &amp; Molecular Biology</i> (2019)</a>.')
SUMMARY = ("Permeabilized nuclei receive a three-base first barcode through indexed "
           "RNA reverse transcription or chromatin tagmentation. Three splinted "
           "ligation rounds add BC2–BC4 and a UMI; TdT/anchor amplification is then "
           "split by restriction digest into RNA and ATAC libraries.")
CAVEAT = ("Native variable TdT tail/Anchor extension and first-cycle tagmented/PCR geometry remain unassigned; complete endpoint assembly is a conditional reference. The source P7XX primer ends in C (one base shorter than full canonical Read2); the retained final R04 terminalT is an extension-derived endpoint boundary, not a silently added primer base. Round-4 oligos contain a 0–3-nt phase block. Drawings use a published "
          "zero-phase member so every downstream primer site remains exact; the phase "
          "position is stated in the oligo list.")


def _seg(name, bases, tag=None, *, placeholder=False, molecular_feature=None):
    return Segment(name, bases, tag, placeholder=placeholder,
                   feature=molecular_feature)


BC1 = feature("paired_seq_cell_bc1", "cell_barcode", "whitelist",
              group="paired_seq_cell_barcode", part="round 1",
              whitelist="Zhu et al. Supplementary Table 2, BC#01")
BC2 = feature("paired_seq_cell_bc2", "cell_barcode", "whitelist",
              group="paired_seq_cell_barcode", part="round 2",
              whitelist="Zhu et al. Supplementary Table 2, BC#02")
BC3 = feature("paired_seq_cell_bc3", "cell_barcode", "whitelist",
              group="paired_seq_cell_barcode", part="round 3",
              whitelist="Zhu et al. Supplementary Table 2, BC#03")
BC4 = feature("paired_seq_cell_bc4", "cell_barcode", "whitelist",
              group="paired_seq_cell_barcode", part="round 4",
              whitelist="Zhu et al. Supplementary Table 2, BC#04")
UMI_FEATURE = feature("paired_seq_umi", "umi", "random")
I5_FEATURE = feature("paired_seq_i5", "sample_index", "unknown")
I7_FEATURE = feature("paired_seq_i7", "sample_index", "unknown")

BC1_COMMON_ARM = "AGGCCAGAGCATTCG"
DNA_BC1_PREFIX = BC1_COMMON_ARM + "AG"
RNA_BC1_PREFIX = BC1_COMMON_ARM + "TC"
NOTI = "GCGGCCGC"
SBFI = "CCTGCAGG"
RT_DT = "T" * 16 + "VN"
MENT_BOTTOM = nx.ME_RC

@dataclass(frozen=True)
class PairedBarcodeRound:
    number:int
    substrate_arm:str
    barcode_5_arm:str
    barcode_3_arm:str
    linker:str
    blocking_oligo:str
    identity:object
    def product_segments(self):
        return [_seg(f"round {self.number} outgoing arm",self.barcode_5_arm),_seg(f"round {self.number} barcode","N"*7,"cbc",placeholder=True,molecular_feature=self.identity),_seg(f"round {self.number} ligating arm",self.barcode_3_arm)]
    def adapter_scene(self):
        if self.linker!=revcomp(self.substrate_arm)+revcomp(self.barcode_3_arm) or self.blocking_oligo!=revcomp(self.linker):raise ValueError("Published splint/blocker relation failed")
        sc=Scene();sc.strand("linker",[_seg("previous splint",revcomp(self.substrate_arm)),_seg("incoming splint",revcomp(self.barcode_3_arm))],label=f"LinkerR0{self.number}")
        sc.anneal("previous",[_seg("previous5′arm",self.substrate_arm)],to="linker",pair=("previous5′arm","previous splint"),label="preceding molecule",mod5="/5Phos/")
        sc.anneal("barcode",self.product_segments(),to="linker",pair=(f"round {self.number} ligating arm","incoming splint"),label=f"BC0{self.number}",mod5="/5Phos/")
        sc.ligation_end("previous","previous5′arm",kind="splinted",polarity="5-prime",owner="substrate",text="source5′P at nick")
        sc.ligation_end("barcode",f"round {self.number} ligating arm",kind="splinted",polarity="3-prime",owner="adapter",text="source3′OH aligned by linker")
        return sc
ROUND2=PairedBarcodeRound(2,BC1_COMMON_ARM,"GTGCGAACTCAGACC","ATCCACGTGCTTGAG","CGAATGCTCTGGCCTCTCAAGCACGTGGAT","ATCCACGTGCTTGAGAGGCCAGAGCATTCG",BC2)
ROUND3=PairedBarcodeRound(3,"GTGCGAACTCAGACC","CATCGGCGTACGACT","GGATTCGAGGAGCGT","GGTCTGAGTTCGCACACGCTCCTCGAATCC","GGATTCGAGGAGCGTGTGCGAACTCAGACC",BC3)

R4_READ2_ARM = "CAGACGTGTGCTCTTCCGATCT"
R4_LIGATING_ARM = "GTGGCCGATGTTTCG"
LINKER_R04 = "AGTCGTACGCCGATGCGAAACATCGGCCAC"
if LINKER_R04 != revcomp(ROUND3.barcode_5_arm) + revcomp(R4_LIGATING_ARM):
    raise ValueError("Paired-seq round-4 linker does not bridge rounds 3 and 4")

ANCHOR_HANDLE = "AAGCAGTGGTATCAACGCAGAGTGAAT"
ANCHOR_G = "G" * 9
ANCHOR_WRITTEN = ANCHOR_HANDLE + ANCHOR_G + "H"
PA_F = R4_READ2_ARM
PA_R = "AAGCAGTGGTATCAACGCAGAGT"
N5XX = il.P5 + "N" * 8 + nx.S5
P7XX = il.P7 + "N" * 6 + il.TRUSEQ_READ2[:-1]
# R04's 22-nt Read-2 arm is the 3' end of the P7XX tail; P7XX adds the remainder.
if not revcomp(il.TRUSEQ_READ2).startswith(revcomp(R4_READ2_ARM)):
    raise ValueError("P7XX tail does not overlap the round-4 Read-2 arm")
P7XX_EXTENSION = revcomp(il.TRUSEQ_READ2)[len(R4_READ2_ARM):]


def _r4_segments() -> list[Segment]:
    return [
        _seg("round 4 Read-2 arm", R4_READ2_ARM, "r2"),
        _seg("round 4 UMI", "U" * 10, "umi", placeholder=True,
             molecular_feature=UMI_FEATURE),
        _seg("round 4 barcode", "D" * 7, "cbc", placeholder=True,
             molecular_feature=BC4),
        _seg("round 4 ligating arm", R4_LIGATING_ARM),
    ]


def _base_segments(kind: str) -> list[Segment]:
    if kind == "RNA":
        return [
            _seg("BC1 common arm", RNA_BC1_PREFIX),
            _seg("round 1 barcode", "A" * 3, "cbc", placeholder=True,
                 molecular_feature=BC1),
            _seg("SbfI site", SBFI),
            _seg("oligo-dT/VN", RT_DT, placeholder=True),
            _seg("first-strand cDNA", "X" * 34, placeholder=True),
        ]
    if kind == "ATAC":
        return [
            _seg("BC1 common arm", DNA_BC1_PREFIX),
            _seg("round 1 barcode", "A" * 3, "cbc", placeholder=True,
                 molecular_feature=BC1),
            _seg("NotI site", NOTI),
            _seg("mosaic end", nx.ME, "me"),
            _seg("accessible genomic DNA", "X" * 40, placeholder=True),
        ]
    raise ValueError(f"unknown Paired-seq modality {kind!r}")


def _prefix(through_round: int) -> list[Segment]:
    rounds = []
    if through_round >= 2:
        rounds.append(ROUND2.product_segments())
    if through_round >= 3:
        rounds.append(ROUND3.product_segments())
    if through_round >= 4:
        rounds.append(_r4_segments())
    return [s for group in reversed(rounds) for s in group]


def _indexed_construct(kind: str, through_round: int) -> Construct:
    return Construct([*_prefix(through_round), *_base_segments(kind)],
                     name=f"Paired-seq {kind} through barcode round {through_round}")


def _indexed_scene(kind: str, through_round: int) -> Scene:
    con = _indexed_construct(kind, through_round)
    sc = Scene()
    sc.strand(kind.lower(), list(con), label=f"{kind} indexed strand",
              mod5="/5Phos/" if through_round < 4 else "")
    previous = "BC1 common arm"
    if through_round >= 2:
        sc.junction(kind.lower(), "round 2 ligating arm", previous,
                    "round 2 ligation")
        previous = "round 2 outgoing arm"
    if through_round >= 3:
        sc.junction(kind.lower(), "round 3 ligating arm", previous,
                    "round 3 ligation")
        previous = "round 3 outgoing arm"
    if through_round >= 4:
        sc.junction(kind.lower(), "round 4 ligating arm", previous,
                    "round 4 ligation")
    sc.labels(kind.lower())
    return sc


def _both_modalities(through_round: int) -> list[Row]:
    return [*_indexed_scene("RNA", through_round).rows(),
            *_indexed_scene("ATAC", through_round).rows()]


def _round4_geometry() -> AdapterLigationGeometry:
    linker = [
        _seg("round 4 previous-arm splint", revcomp(ROUND3.barcode_5_arm)),
        _seg("round 4 barcode splint", revcomp(R4_LIGATING_ARM)),
    ]
    sc = Scene()
    sc.strand("round 4 linker", linker, label="Linker-R04")
    sc.anneal("round 3 tagged molecule",
              [_seg("round 3 terminal 5-prime arm", ROUND3.barcode_5_arm)],
              to="round 4 linker",
              pair=("round 3 terminal 5-prime arm",
                    "round 4 previous-arm splint"),
              label="round 3 tagged molecule", mod5="/5Phos/")
    sc.anneal("round 4 barcode oligo", _r4_segments(), to="round 4 linker",
              pair=("round 4 ligating arm", "round 4 barcode splint"),
              label="BC#04 oligo",
              unpaired=("round 4 Read-2 arm", "round 4 UMI",
                        "round 4 barcode"))
    sc.ligation_end("round 3 tagged molecule",
                    "round 3 terminal 5-prime arm",
                    kind="splinted", polarity="5-prime", owner="substrate",
                    text="5′-phosphate at the ligation nick")
    sc.ligation_end("round 4 barcode oligo", "round 4 ligating arm",
                    kind="splinted", polarity="3-prime", owner="adapter",
                    text="3′-OH at the ligation nick")
    sc.labels("round 4 linker")
    sc.labels("round 3 tagged molecule")
    sc.labels("round 4 barcode oligo")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Supplementary Tables 1–2 define both linker halves, the 10-nt UMI, BC4 and the variable phase block.")


def _geometry(round_no: int) -> AdapterLigationGeometry:
    round_ = {2: ROUND2, 3: ROUND3}[round_no]
    return AdapterLigationGeometry.known(
        "splinted", round_.adapter_scene(),
        statement=f"Supplementary Tables 1–2 define both arms of Linker-R0{round_no} and the complete barcode oligo.")


def _tailed_scene(kind):
    return molecule(kind+" actual tailed ssDNA",tuple([*list(_indexed_construct(kind,4)),role("TdT poly(C) tract; variable length")]))

def _anchor_scene(kind):
    sc=_tailed_scene(kind)
    sc.strand("anchor",[_seg("source Anchor",ANCHOR_WRITTEN,placeholder=True)],label="independent added Anchor; native annealing register/extension extent unassigned")
    return sc

def _final_core(kind: str) -> list[Segment]:
    if kind == "RNA":
        branch = [_seg("poly(A) copied from oligo-dT/VN", "NB" + "A" * 16,
                       placeholder=True),
                  _seg("SbfI site", SBFI),
                  _seg("round 1 barcode", "A" * 3, "cbc", placeholder=True,
                       molecular_feature=BC1),
                  _seg("RNA BC1 exit arm", "GACGAATGCTCTGGCCT")]
        payload = _seg("cDNA", "X" * 42, placeholder=True)
    else:
        branch = [_seg("right mosaic end reverse complement", nx.ME_RC, "me"),
                  _seg("NotI site", NOTI),
                  _seg("round 1 barcode", "A" * 3, "cbc", placeholder=True,
                       molecular_feature=BC1),
                  _seg("ATAC BC1 exit arm", "CTCGAATGCTCTGGCCT")]
        payload = _seg("accessible genomic DNA", "X" * 42, placeholder=True)
    return [
        _seg("P5", il.P5, "p5"),
        _seg("i5", "J" * 8, "cbc", placeholder=True,
             molecular_feature=I5_FEATURE),
        _seg("S5", nx.S5, "s5"),
        _seg("left mosaic end", nx.ME, "me"),
        payload,
        *branch,
        _seg("round 2 incoming arm", revcomp(ROUND2.barcode_3_arm)),
        _seg("round 2 barcode", "B" * 7, "cbc", placeholder=True,
             molecular_feature=BC2),
        _seg("round 2 outgoing arm", revcomp(ROUND2.barcode_5_arm)),
        _seg("round 3 incoming arm", revcomp(ROUND3.barcode_3_arm)),
        _seg("round 3 barcode", "C" * 7, "cbc", placeholder=True,
             molecular_feature=BC3),
        _seg("round 3 outgoing arm", revcomp(ROUND3.barcode_5_arm)),
        _seg("round 4 incoming arm", revcomp(R4_LIGATING_ARM)),
        _seg("round 4 barcode", "D" * 7, "cbc", placeholder=True,
             molecular_feature=BC4),
        _seg("round 4 UMI", "U" * 10, "umi", placeholder=True,
             molecular_feature=UMI_FEATURE),
        _seg("round 4 Read-2 arm", revcomp(R4_READ2_ARM), "r2"),
        _seg("P7XX adapter extension", P7XX_EXTENSION, "r2"),
        _seg("i7 reverse complement", "I" * 6, "cbc", placeholder=True,
             molecular_feature=I7_FEATURE),
        _seg("P7 reverse complement", il.P7_RC, "p7"),
    ]


PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def _library(kind: str) -> Construct:
    lib = Construct([replace(part,inferred=True) for part in _final_core(kind)], name=f"Paired-seq {kind} library")
    problems = sp.verify(lib, PRIMERS,
                         required_roles=tuple(p.role for p in PRIMERS))
    if problems:
        raise ValueError(f"invalid Paired-seq {kind} library: " + "; ".join(problems))
    return lib


RNA = _library("RNA")
ATAC = _library("ATAC")
FINAL_LIBRARIES = (
    ("RNA library", RNA, PRIMERS,
     "INFERRED endpoint — zero-phase representative. Read 1 enters cDNA; Read 2 begins with UMI and the four-part cell barcode chain.",
     "The published mixed primer set is derived from the completed library."),
    ("Chromatin-accessibility library", ATAC, PRIMERS,
     "INFERRED endpoint — zero-phase representative. Read 1 enters accessible DNA; Read 2 begins with UMI and the matched four-part cell barcode chain.",
     "The same published mixed primer set is derived from the completed library."),
)
READ_LENGTHS = {
    "RNA library": {"Read 1": 53, "Index 1 (i7)": 7, "Read 2": 130},
    "Chromatin-accessibility library": {
        "Read 1": 53, "Index 1 (i7)": 7, "Read 2": 130},
}
READ_USES = {
    "RNA library": {"Read 1": "include", "Index 1 (i7)": "demultiplex",
                    "Read 2": "extract"},
    "Chromatin-accessibility library": {
        "Read 1": "include", "Index 1 (i7)": "demultiplex",
        "Read 2": "extract"},
}


def oligos():
    return (
        f"<b>Barcoded Tn5 top</b> — 5′-/5Phos/{DNA_BC1_PREFIX}[3-nt BC1]{NOTI}{nx.ME}-3′",
        f"<b>AdaptorA (printed; exact Tn5-P5 load unassigned)</b> — 5′-{nx.READ1_PRIMER}-3′",
        f"<b>pMENTs</b> — 5′-/5Phos/{MENT_BOTTOM}/ddC/-3′",
        f"<b>RT primer</b> — 5′-/5Phos/{RNA_BC1_PREFIX}[3-nt BC1]{SBFI}{RT_DT}-3′",
        f"<b>BC2</b> — 5′-/5Phos/{ROUND2.barcode_5_arm}[7-nt BC2]{ROUND2.barcode_3_arm}-3′",
        f"<b>Linker-R02</b> — 5′-{ROUND2.linker}-3′",
        f"<b>Blocker-R02</b> — 5′-{ROUND2.blocking_oligo}-3′",
        f"<b>BC3</b> — 5′-/5Phos/{ROUND3.barcode_5_arm}[7-nt BC3]{ROUND3.barcode_3_arm}-3′",
        f"<b>Linker-R03</b> — 5′-{ROUND3.linker}-3′",
        f"<b>Blocker-R03</b> — 5′-{ROUND3.blocking_oligo}-3′",
        f"<b>BC4</b> — 5′-{R4_READ2_ARM}[10-nt UMI][7-nt BC4][0–3-nt phase]{R4_LIGATING_ARM}-3′",
        f"<b>Linker-R04</b> — 5′-{LINKER_R04}-3′",
        f"<b>Terminator</b> — 5′-{revcomp(LINKER_R04)}-3′",
        f"<b>Anchor</b> — 5′-{ANCHOR_WRITTEN}-3′",
        f"<b>PA-F / PA-R</b> — 5′-{PA_F}-3′ / 5′-{PA_R}-3′",
        f"<b>N5XX</b> — 5′-{il.P5}[8-nt i5]{nx.S5}-3′",
        f"<b>P7XX</b> — 5′-{il.P7}[6-nt i7]{il.TRUSEQ_READ2[:-1]}-3′",
    )


def _initial_rows() -> list[Row]:
    chromatin = Scene.duplex([_seg("fixed chromatin", "X" * 40,
                                  placeholder=True)], label="fixed chromatin")
    mrna = Scene()
    mrna.strand("mRNA", [_seg("transcript", "X" * 28, placeholder=True),
                         _seg("poly(A)", "A" * 16)], label="poly(A) RNA")
    return [*chromatin.rows(), *mrna.rows()]


def _single_strand(name: str, bases: str) -> Scene:
    sc = Scene()
    sc.strand(name, [_seg(name, bases)], label=name)
    return sc


def mixed_state(label,*,rounds=4,tailed=False,extended=False):
    sc=Scene()
    for kind in ("RNA","ATAC"):
        parts=list(_indexed_construct(kind,rounds))
        if tailed:parts.append(role("TdT poly(C) tract; variable length"))
        if extended:
            # Retain the source-known ancestral barcode chain without assigning
            # a tail length or a complete newly extended strand.
            parts += [role("variable TdT tail and native extension-end assembly unassigned",inferred=True)]
        sc.strand(kind,parts,label=label+" independent "+kind+" identifier-bearing strand class")
        if extended:sc.strand(kind+" newcopies",[role(kind+" newly extended/copied complementary strands; complete extent unassigned",inferred=True)],label="independent source-generated copy population")
    return sc

def cleanup(g,key,template,recipe,library=False):
    rule=g.input(key+"_rule",recipe,MaterialKind.SELECTION)
    out=g.product(key+"_clean","purified "+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(key+"_waste","source column/bead/size excluded fractions",MaterialKind.WASTE)
    g.react(key,"Purify source product",StepKind.CLEANUP,((template,InputRole.TEMPLATE),(rule,InputRole.SELECTION)),(out,waste),changes_molecule=False);g.discard(waste,"Source physical bead/column/wash transfer; molecular cleavage alone is not removal")
    return out

def retagmented_state(prefix,prior):
    sc=Scene()
    sc.strand("tagged target",[*list(_indexed_construct("RNA" if prefix=="rna" else "ATAC",4))[:-1],role(prefix+" fragment interval after retagmentation; extent unassigned",inferred=True),role("new P5-associated end class; native strand/gap assembly unassigned",inferred=True)],label="INFERRED barcode-retaining native pre-index target class")
    sc.strand("carried cleavage products",[role("carried opposite-modality cleavage fragments and other tagmented classes; recovery unassigned")],label="independent carried molecular population")
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nuclei=g.input("nuclei","permeabilized nuclei; independent chromatin and RNA",MaterialKind.BIOLOGICAL,state=_initial_rows(),primary=True)
    trans=Scene();trans.strand("top",_base_segments("ATAC")[:-1],label="sourceBC1 transferred Tn5strand",mod5="/5Phos/")
    # pMENTs ddC is a physical terminal nucleotide on the nontransferred oligo.
    trans.anneal("bottom",[Segment("pMENTs",MENT_BOTTOM),Segment("terminalddC","C",note="3′dideoxyC; not copied into transferred strand")],to="top",pair=("pMENTs","mosaic end"),unpaired=("terminalddC",),label="sourcepMENTs nontransferred strand",mod5="/5Phos/",mod3="ddC")
    adapter=g.input("bc1_transposome","source barcodedTn5 top/pMENTs; actual annealed geometry bounded",MaterialKind.ADAPTER,state=trans)
    tn5=g.input("tn5","loadedTn5;37°C30min,EDTAstop",MaterialKind.CATALYST)
    ts=Scene();ts.strand("ATAC",list(_indexed_construct("ATAC",1)),label="transferred BC1DNA strand; native duplex/gaps unassigned");ts.strand("RNA",[role("cellular RNA")],label="independent unchangedRNA")
    tagged=g.product("tagged","nuclei with tagged DNA and RNA",MaterialKind.MOLECULAR,state=ts)
    other=g.product("tn5_other","other transposition end classes",MaterialKind.MOLECULAR)
    g.react("tagment_chromatin","Transpose chromatin with matchedBC1",StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Native end-class amplification compatibility does not prove physical discard; recovery through nuclei pellet unassigned."))
    tagged=cleanup(g,"tagment_recovery",tagged,"1000g10min4°C; source removes supernatant and resuspends nuclei")
    rp=g.input("rt_primer","source matchedBC1 oligo(dT16)VN",MaterialKind.OLIGO,state=Construct(_base_segments("RNA")[:-1]))
    rt=g.input("rt","MaximaEP0743;source50°C/rampcycling",MaterialKind.CATALYST)
    copied=g.product("copied","actual BC1ss-cDNA beside taggedDNA",MaterialKind.MOLECULAR,state=mixed_state("RT",rounds=1))
    g.react("reverse_transcribe","Reverse-transcribe RNA with matchedBC1 primer",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(rp,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
    current=cleanup(g,"rt_recovery",copied,"sourcepool/Triton/1000g10min4°C and supernatant removal")
    ligase=g.input("barcode_ligase","T4DNAligaseM0202L;37°C30min",MaterialKind.CATALYST)
    for number,geometry in ((2,_geometry(2)),(3,_geometry(3)),(4,_round4_geometry())):
        # Reagent-only inputs do not include the preceding substrate strand.
        reagent=Scene()
        if number<4:
            definition=ROUND2 if number==2 else ROUND3
            reagent.strand("barcode",definition.product_segments(),label="5′phosphorylated barcode",mod5="/5Phos/")
            reagent.strand("splint",[Segment("linker",definition.linker)],label="independent splint")
            blocker=definition.blocking_oligo
        else:
            reagent.strand("barcode",_r4_segments(),label="sourceR04barcode")
            reagent.strand("splint",[Segment("linker",LINKER_R04)],label="independent LinkerR04")
            blocker=revcomp(LINKER_R04)
        barcode=g.input(f"barcode_{number}",f"source barcode/linker round{number}",MaterialKind.OLIGO,state=reagent)
        stop=g.input(f"blocker_{number}","source blocker/terminator added after30min;R04includesEDTA",MaterialKind.OLIGO,state=blocker)
        joined=g.product(f"round_{number}",f"independent RNA/DNA afterBC{number}",MaterialKind.MOLECULAR,state=mixed_state(f"BC{number}",rounds=number))
        g.react(f"ligate_{number}",f"Ligate BC{number} then add blocker/terminator",StepKind.LIGATION,((current,InputRole.TEMPLATE),(barcode,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST),(stop,InputRole.STRUCTURAL)),(joined,),geometry=geometry)
        current=joined
    rule=g.input("sublibrary_rule","count nuclei; source10k–100k nuclei/sublibrary, optimal~25k",MaterialKind.CONDITION)
    chosen=g.product("chosen_sublibrary","chosen mixed nuclei sublibrary",MaterialKind.MOLECULAR,state=current.state)
    rest=g.product("other_sublibraries","other nuclei sublibraries",MaterialKind.MOLECULAR,state=current.state)
    g.react("make_sublibraries","Divide barcoded nuclei into sublibraries",StepKind.SPLIT,((current,InputRole.TEMPLATE),(rule,InputRole.CONDITION)),(chosen,rest),changes_molecule=False);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,"Other sublibrary storage/use unassigned in this representation."))
    lysis=g.input("lysis","NaCl/SDS/PK55°C2h",MaterialKind.CONDITION)
    released=g.product("released","released source DNA/cDNA mixed population",MaterialKind.MOLECULAR,state=chosen.state)
    g.react("lyse","Lyse selected nuclei",StepKind.DEPROTEINIZATION,((chosen,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(released,))
    released=cleanup(g,"lysis_cleanup",released,"1×SPRI;12.5µlwater")
    den=g.input("denature_condition","95°C5min,ice5min",MaterialKind.CONDITION)
    ss=g.product("denatured","denatured DNA/cDNA strand populations",MaterialKind.MOLECULAR,state=released.state)
    g.react("denature","Denature before TdT",StepKind.THERMAL_TREATMENT,((released,InputRole.TEMPLATE),(den,InputRole.CONDITION)),(ss,))
    tdt=g.input("tdt","TdT+dCTP;37°C30min,75°C20min",MaterialKind.CATALYST)
    tailed=g.product("tailed","actual variable poly(C)-tailed strands",MaterialKind.MOLECULAR,state=mixed_state("TdT",tailed=True))
    g.react("tail","Add variable poly(C) tract",StepKind.TAILING,((ss,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
    anchor=g.input("anchor","source Anchor9G/H primer",MaterialKind.OLIGO,state=ANCHOR_WRITTEN)
    extpol=g.input("extension_pol","KAPAHiFiHS;15source cycles with two47/68°C intervals percycle",MaterialKind.CATALYST)
    extended=g.product("extended","Anchor-extension target classes; native length/first-copy geometry unassigned",MaterialKind.MOLECULAR,state=mixed_state("Anchor extension",extended=True))
    g.react("extend_anchor","Extend/amplify source Anchor15cycles",StepKind.EXTENSION,((tailed,InputRole.TEMPLATE),(anchor,InputRole.STRUCTURAL),(extpol,InputRole.CATALYST)),(extended,))
    pf=g.input("pa_f","source PA-F",MaterialKind.OLIGO,state=PA_F)
    pr=g.input("pa_r","source PA-R",MaterialKind.OLIGO,state=PA_R)
    preamp=g.product("preamp","preamplified mixed RNA/DNA products",MaterialKind.MOLECULAR,state=mixed_state("PA preamp",extended=True))
    g.react("preamp_pcr","Preamplify PA-F/PA-R10cycles",StepKind.PCR,((extended,InputRole.TEMPLATE),(pf,InputRole.STRUCTURAL),(pr,InputRole.STRUCTURAL),(extpol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(extended.state,reason="Actual Anchor-extension classes enter PCR. Tail length and native complete extension/handle assembly unassigned; fixed9C final/exact template not substituted."))
    preamp=cleanup(g,"preamp_cleanup",preamp,"SPRI double-size10µl+32.5µl;34µlelution")
    qc=g.product("preamp_qc","1µl preamplified product for quantification",MaterialKind.MOLECULAR,state=preamp.state)
    ra=g.product("rna_aliquot","RNA-preparation aliquot; both modalities retained",MaterialKind.MOLECULAR,state=preamp.state)
    aa=g.product("atac_aliquot","ATAC-preparation aliquot; both modalities retained",MaterialKind.MOLECULAR,state=preamp.state)
    g.react("split","Take1µlQC; divide remaining33µl into two mixed aliquots",StepKind.SPLIT,((preamp,InputRole.TEMPLATE),),(qc,ra,aa),changes_molecule=False);g.mark_terminal(Terminal(qc.id,Disposition.UNRESOLVED,"Source quantification aliquot; subsequent fate unassigned."))
    for prefix,template,final,enzyme,heading in (("rna",ra,RNA,"NotI-HF",FINAL_LIBRARIES[0][0]),("atac",aa,ATAC,"SbfI-HF",FINAL_LIBRARIES[1][0])):
        re=g.input(prefix+"_restriction",enzyme+";37°C60min",MaterialKind.CATALYST)
        state=Scene();state.strand("target",[role(prefix+" target molecule not cut at opposite-modality site; native ends unassigned")],label="surviving target architecture")
        state.strand("cut_other",[role("restriction-cleaved opposite-modality fragments; not physically purified by nuclease")],label="independent cleavage-product population")
        cut=g.product(prefix+"_cut","digested mixed "+prefix+" preparation",MaterialKind.MOLECULAR,state=state)
        g.react(prefix+"_digest","Digest opposite modality with "+enzyme,StepKind.DIGESTION,((template,InputRole.TEMPLATE),(re,InputRole.CATALYST)),(cut,))
        cut=cleanup(g,prefix+"_restriction_cleanup",cut,"sourceQIAquickPCRpurification;30µl0.1×EB")
        adapter=g.input(prefix+"_p5_adapter","source Tn5-P5; exact loadedadapter composition unassigned (Table1 prints AdaptorA without assigning this load)",MaterialKind.ADAPTER)
        tn=g.input(prefix+"_tn5","Tn5-P5;37°C30min",MaterialKind.CATALYST)
        tagged=g.product(prefix+"_tagged","actual pre-index tagmented mixed "+prefix+" preparation",MaterialKind.MOLECULAR,state=retagmented_state(prefix,state))
        other=g.product(prefix+"_other_ends","other Tn5 end classes",MaterialKind.MOLECULAR)
        g.react(prefix+"_tagment","Retagment source preparation withTn5-P5",StepKind.TAGMENTATION,((cut,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn,InputRole.CATALYST)),(tagged,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Nonamplifiable end classes are not asserted physically removed; source cleanup recovery unassigned."))
        tagged=cleanup(g,prefix+"_tagment_cleanup",tagged,"QIAquickPCRpurification;30µl0.1×EB")
        f=g.input(prefix+"_n5","source N5XX8N index primer",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),Segment("i5","N"*8,placeholder=True,feature=I5_FEATURE),Segment("S5",nx.S5)]))
        r=g.input(prefix+"_p7","source P7XX6N index primer; printed32nt Read2arm",MaterialKind.OLIGO,state=Construct([Segment("P7",il.P7),Segment("i7","N"*6,placeholder=True,feature=I7_FEATURE),Segment("source Read2arm",il.TRUSEQ_READ2[:-1])]))
        q5=g.input(prefix+"_q5","NEBQ5;72°C5min, source Step2 additional10–15repeats to10nM",MaterialKind.CATALYST)
        amp=g.product(prefix+"_amp","conditional "+prefix+" endpoint component",MaterialKind.MOLECULAR,state=final,inferred=True)
        g.react(prefix+"_index_pcr","INFERRED native assembly — complete indexed target",StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(tagged.state,reason="Actual tagged mixed input retained. Native gap filling/first-cycle copying and complete primer-to-product junctions unassigned; endpoint is conditional reference."))
        lib=cleanup(g,prefix+"_final",amp,"0.85×SPRI",library=True)
        g.sequence(lib,SequencingPlan(library=final,primers=PRIMERS,acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths=READ_LENGTHS[heading],read_uses=READ_USES[heading]))
    g.annotate(Highlight("reverse_transcribe","RNA-template fate unassigned","RT adds cDNA; subsequent nativeRNA recovery/degradation is not equated with physical removal."))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    return g.require_publishable()

SHOW_READ_CONTENT=False
READ_SPECS={name:(sp.read_spec("Read 1",53,"biological insert","include"),sp.read_spec("Index 1 (i7)",7,"6nt printed index plus adjacent primer base; seventh-cycle interpretation derived","demultiplex"),sp.read_spec("Read 2",130,"UMI and four-coordinate cell barcode chain","extract"),sp.read_spec("Index 2 (i5)","unassigned","physical8nt primer index; acquisition not assigned by stated53+7+130 run","undetermined")) for name,*_ in FINAL_LIBRARIES}
READ_SPECS_NOTE="Positive published53+7+130 run retained. Physical8nt i5 primer window does not establish index2 acquisition."
REVIEWED_TARGETS=("nuclei","bc1_transposome","tn5","tagment_chromatin","tagment_recovery_rule","tagment_recovery","rt_primer","rt","reverse_transcribe","rt_recovery_rule","rt_recovery","barcode_ligase","barcode_2","blocker_2","ligate_2","barcode_3","blocker_3","ligate_3","barcode_4","blocker_4","ligate_4","sublibrary_rule","make_sublibraries","lysis","lyse","lysis_cleanup_rule","lysis_cleanup","denature_condition","denature","tdt","tail","anchor","extension_pol","extend_anchor","pa_f","pa_r","preamp_pcr","preamp_cleanup_rule","preamp_cleanup","split","rna_restriction","rna_digest","rna_restriction_cleanup_rule","rna_restriction_cleanup","rna_p5_adapter","rna_tn5","rna_tagment","rna_tagment_cleanup_rule","rna_tagment_cleanup","rna_n5","rna_p7","rna_q5","rna_index_pcr","rna_final_rule","rna_final","rna_final_clean","atac_restriction","atac_digest","atac_restriction_cleanup_rule","atac_restriction_cleanup","atac_p5_adapter","atac_tn5","atac_tagment","atac_tagment_cleanup_rule","atac_tagment_cleanup","atac_n5","atac_p7","atac_q5","atac_index_pcr","atac_final_rule","atac_final","atac_final_clean")
