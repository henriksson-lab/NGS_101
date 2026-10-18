"""Source-backed PairedTag DNA/RNA chemistry.

Representative plate members are from Zhu et al. Supplementary Table 1.  The
round-3 plate contains 0--3 phase bases; molecular states use its zero-phase A1
member so the barcode and primer geometry remain exact.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
from dataclasses import dataclass,replace
from protocol_model import molecule,role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind, SequencingPlan, Disposition, Terminal, Highlight)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


TITLE = "PairedTag — joint histone-mark and transcriptome profiling"
NOTES = "01_pairedtag.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41592-021-01060-3">'
          'Zhu et al., <i>Nature Methods</i> (2021)</a>.')
SUMMARY = ("A BC1-loaded pA–Tn5 marks antibody-bound chromatin while matched BC1 "
           "oligo(dT) and random-hexamer primers copy nuclear RNA. Two exact splinted "
           "ligations add BC2, BC3 and a UMI before modality-specific library finishing.")
CAVEAT = ("Round-3 barcode oligos contain a 0–3-nt phase block. Drawings use the "
          "published zero-phase A1 member. The four-base FokI product overhang varies "
          "with the cut molecule; the adapter's source-printed NNDC cohesive end is shown.")

# Supplementary Table 1 common oligos.
PMENTS = nx.ME_RC
ADAPTOR_A = nx.ADAPTOR_S5
LINKER_R02 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
BLOCKER_R02 = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
LINKER_R03 = "GGTCTGAGTTCGCACCGAAACATCGGCCAC"
QUENCHER_R03 = "GTGGCCGATGTTTCGGTGCGAACTCAGACC"

BC1_COMMON = "AGGCCAGAGCATTCG"
DNA_DISCRIMINATOR = "A"
RNA_DISCRIMINATOR = "T"
NOTI = "GCGGCCGC"
SBFI = "CCTGCAGG"
RNA_DT = "T" * 16 + "VN"
RNA_N6 = "N" * 6

@dataclass(frozen=True)
class TagRound2:
    substrate_arm:str
    barcode_5_arm:str
    barcode_3_arm:str
    linker:str
    feature:object
    def product_segments(self):
        return [seg("round 2 outgoing arm",self.barcode_5_arm),seg("round 2 barcode","N"*7,"cbc",placeholder=True,feature=self.feature),seg("round 2 ligating arm",self.barcode_3_arm)]
    def adapter_scene(self):
        if self.linker!=revcomp(self.substrate_arm)+revcomp(self.barcode_3_arm):raise ValueError("Published PairedTag splint arms do not match")
        sc=Scene();sc.strand("linker",[seg("previous splint",revcomp(self.substrate_arm)),seg("incoming splint",revcomp(self.barcode_3_arm))],label="LinkerR02")
        sc.anneal("preceding",[seg("previous5′arm",self.substrate_arm)],to="linker",pair=("previous5′arm","previous splint"),label="preceding molecule",mod5="/5Phos/")
        sc.anneal("barcode",self.product_segments(),to="linker",pair=("round 2 ligating arm","incoming splint"),label="sourceBC2",mod5="/5Phos/")
        sc.ligation_end("preceding","previous5′arm",kind="splinted",polarity="5-prime",owner="substrate",text="source5′P")
        sc.ligation_end("barcode","round 2 ligating arm",kind="splinted",polarity="3-prime",owner="adapter",text="source3′OH")
        return sc
ROUND2=TagRound2(BC1_COMMON,"GTGCGAACTCAGACC","ATCCACGTGCTTGAG",LINKER_R02,feature("pairedtag_cell_bc2","cell_barcode","whitelist",group="pairedtag_cell_barcode",part="round2",whitelist="publishedTable1BCplate02"))

R3_READ2_ARM = "CAGACGTGTGCTCTTCCGATCT"
R3_LIGATING_ARM = "GTGGCCGATGTTTCG"
R3_UMI_LENGTH = 10
R3_BARCODE_LENGTH = 7
if LINKER_R03 != revcomp(ROUND2.barcode_5_arm) + revcomp(R3_LIGATING_ARM):
    raise ValueError("Linker-R03 must bridge the round-2 and round-3 oligos")
if BLOCKER_R02 != revcomp(LINKER_R02):
    raise ValueError("Blocker-R02 must be the Linker-R02 reverse complement")
if QUENCHER_R03 != revcomp(LINKER_R03):
    raise ValueError("Quencher-R03 must be the Linker-R03 reverse complement")

ANCHOR_HANDLE = "AAGCAGTGGTATCAACGCAGAGTGAAGGATGT"
ANCHOR_G = "G" * 9
ANCHOR_WRITTEN = ANCHOR_HANDLE + ANCHOR_G + "H"
P5_FOKI = il.TRUSEQ_READ1
P5C_NNDC_CORE = revcomp(P5_FOKI)[:-1]
P5C_NNDC = "NNDC" + P5C_NNDC_CORE
P5H_FOKI = P5_FOKI + "H"
P5HC_NNDC = "NNDCD" + P5C_NNDC_CORE
PA_F = R3_READ2_ARM
PA_R = "AAGCAGTGGTATCAACGCAGAGT"
N5XX = il.P5 + "N" * 8 + nx.S5
P7XX = il.P7 + "N" * 6 + il.TRUSEQ_READ2[:-1]
P5_UNIVERSAL = il.TRUSEQ_P5_FULL

BC1 = feature("pairedtag_cell_bc1", "cell_barcode", "whitelist",
              group="pairedtag_cell_barcode", part="round 1",
              whitelist="Zhu et al. Supplementary Table 1, BC Plate #01")
BC2 = ROUND2.feature
BC3 = feature("pairedtag_cell_bc3", "cell_barcode", "whitelist",
              group="pairedtag_cell_barcode", part="round 3",
              whitelist="Zhu et al. Supplementary Table 1, BC Plate #03")
UMI = feature("pairedtag_umi", "umi", "random")
I5 = feature("pairedtag_i5", "sample_index", "unknown")
I7 = feature("pairedtag_i7", "sample_index", "unknown")


def _bc1_segments(kind: str) -> list[Segment]:
    if kind == "RNA":
        return [
            seg("BC1 common arm", BC1_COMMON),
            seg("RNA discriminator", RNA_DISCRIMINATOR),
            seg("round-1 barcode", "A" * 4, "cbc", placeholder=True, feature=BC1),
            seg("SbfI site", SBFI),
            seg("oligo-dT16", "T" * 16),
            seg("VN anchor", "VN", placeholder=True),
            seg("first-strand cDNA", "X" * 34, placeholder=True),
        ]
    if kind == "RNA-N6":
        return [
            seg("BC1 common arm", BC1_COMMON),
            seg("RNA discriminator", RNA_DISCRIMINATOR),
            seg("round-1 barcode", "A" * 4, "cbc", placeholder=True, feature=BC1),
            seg("SbfI site", SBFI),
            seg("random hexamer", RNA_N6, placeholder=True),
            seg("first-strand cDNA", "X" * 34, placeholder=True),
        ]
    if kind == "DNA":
        return [
            seg("BC1 common arm", BC1_COMMON),
            seg("DNA discriminator", DNA_DISCRIMINATOR),
            seg("round-1 barcode", "A" * 4, "cbc", placeholder=True, feature=BC1),
            seg("NotI site", NOTI),
            seg("mosaic end", nx.ME, "me"),
            seg("histone-mark DNA", "X" * 40, placeholder=True),
        ]
    raise ValueError(f"unknown PairedTag branch {kind!r}")


def _r3_segments() -> list[Segment]:
    return [
        seg("round-3 Read-2 arm", R3_READ2_ARM, "r2"),
        seg("round-3 UMI", "U" * R3_UMI_LENGTH, "umi", placeholder=True, feature=UMI),
        seg("round-3 barcode", "C" * R3_BARCODE_LENGTH, "cbc",
            placeholder=True, feature=BC3),
        seg("round-3 ligating arm", R3_LIGATING_ARM),
    ]


def _prefix(through_round: int) -> list[Segment]:
    rounds = []
    if through_round >= 2:
        rounds.append(ROUND2.product_segments())
    if through_round >= 3:
        rounds.append(_r3_segments())
    return [s for group in reversed(rounds) for s in group]


def _indexed_construct(kind: str, through_round: int) -> Construct:
    return Construct([*_prefix(through_round), *_bc1_segments(kind)],
                     name=f"PairedTag {kind} through barcode round {through_round}")


def _indexed_scene(kind: str, through_round: int) -> Scene:
    con = _indexed_construct(kind, through_round)
    sc = Scene()
    sc.strand(kind.lower(), list(con), label=f"{kind} indexed strand",
              mod5="/5Phos/" if through_round < 3 else "")
    previous = "BC1 common arm"
    if through_round >= 2:
        sc.junction(kind.lower(), "round 2 ligating arm", previous,
                    "round-2 ligation")
        previous = "round 2 outgoing arm"
    if through_round >= 3:
        sc.junction(kind.lower(), "round-3 ligating arm", previous,
                    "round-3 ligation")
    sc.labels(kind.lower())
    return sc


def _both_modalities(through_round: int) -> list[Row]:
    return [*_indexed_scene("RNA", through_round).rows(),
            *_indexed_scene("DNA", through_round).rows()]


def tn5_loading_scene() -> Scene:
    top = _bc1_segments("DNA")[:-1]
    bottom = [seg("pMENTs mosaic-end complement", PMENTS, "me"),seg("terminal ddC nucleotide","C",note="nontransferred physical3′dideoxyC")]
    sc = Scene()
    sc.strand("barcoded top", top, label="DNA_#01_RE", mod5="/5Phos/")
    sc.anneal("pMENTs", bottom, to="barcoded top",
              pair=("pMENTs mosaic-end complement", "mosaic end"),
              label="pMENTs",unpaired=("terminal ddC nucleotide",), mod5="/5Phos/", mod3="ddC")
    sc.note("pMENTs", "barcoded duplex loaded into pA–Tn5")
    sc.labels("barcoded top")
    return sc


def rt_scene(random: bool = False) -> Scene:
    kind = "RNA-N6" if random else "RNA"
    primer = _bc1_segments(kind)[:-1]
    sc = Scene()
    if random:
        target = [seg("RNA upstream", "N" * 18, placeholder=True),
                  seg("random landing site", "N" * 6, placeholder=True),
                  seg("RNA downstream", "N" * 18, placeholder=True)]
        sc.strand("RNA", target, label="nuclear RNA")
        sc.anneal("RT primer", primer, to="RNA",
                  pair=("random hexamer", "random landing site"),
                  unpaired=("BC1 common arm", "RNA discriminator",
                            "round-1 barcode", "SbfI site"),
                  label="RNA_#01_NRE", mod5="/5Phos/")
    else:
        target = [seg("transcript", "X" * 34, placeholder=True),
                  seg("poly(A)", "A" * 16)]
        sc.strand("RNA", target, label="poly(A) nuclear RNA")
        sc.anneal("RT primer", primer, to="RNA",
                  pair=("oligo-dT16", "poly(A)"),
                  unpaired=("BC1 common arm", "RNA discriminator",
                            "round-1 barcode", "SbfI site", "VN anchor"),
                  label="RNA_#01_RE", mod5="/5Phos/")
    sc.arrow("RT primer", "Maxima H-minus reverse transcriptase")
    sc.labels("RT primer")
    return sc


def _round2_geometry() -> AdapterLigationGeometry:
    return AdapterLigationGeometry.known(
        "splinted", ROUND2.adapter_scene(),
        statement="Linker-R02 pairs the BC1 common arm with the complete BC2 oligo; Blocker-R02 is its exact reverse complement.")


def _round3_geometry() -> AdapterLigationGeometry:
    linker = [
        seg("round-2-arm splint", revcomp(ROUND2.barcode_5_arm)),
        seg("round-3-arm splint", revcomp(R3_LIGATING_ARM)),
    ]
    sc = Scene()
    sc.strand("Linker-R03", linker, label="Linker-R03")
    sc.anneal("BC2-tagged molecule",
              [seg("round-2 terminal 5′ arm", ROUND2.barcode_5_arm)],
              to="Linker-R03",
              pair=("round-2 terminal 5′ arm", "round-2-arm splint"),
              label="BC2-tagged molecule", mod5="/5Phos/")
    sc.anneal("R04 barcode oligo", _r3_segments(), to="Linker-R03",
              pair=("round-3 ligating arm", "round-3-arm splint"),
              unpaired=("round-3 Read-2 arm", "round-3 UMI",
                        "round-3 barcode"), label="R04_#01")
    sc.ligation_end("BC2-tagged molecule", "round-2 terminal 5′ arm",
                    kind="splinted", polarity="5-prime", owner="substrate",
                    text="5′-phosphate at the ligation nick")
    sc.ligation_end("R04 barcode oligo", "round-3 ligating arm",
                    kind="splinted", polarity="3-prime", owner="adapter",
                    text="3′-OH at the ligation nick")
    sc.labels("Linker-R03")
    sc.labels("BC2-tagged molecule")
    sc.labels("R04 barcode oligo")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Linker-R03 aligns the published R04 oligo to BC2; Quencher-R03 is the exact reverse complement and terminates the reaction.")


def _tailed_scene(kind):
    return molecule(kind+" variable TdT tail",tuple([*list(_indexed_construct(kind,3)),role("variable poly(C) tail; length unassigned")]))

def _anchor_scene(kind):
    sc=_tailed_scene(kind);sc.strand("anchor",[seg("sourceAnchor",ANCHOR_WRITTEN,placeholder=True)],label="independent sourceAnchor; annealing register/extension extent unassigned")
    return sc

def p5_foki_geometry() -> AdapterLigationGeometry:
    top = [seg("unpaired 5′ A", P5_FOKI[:1], "r1"),
           seg("P5-FokI paired arm", P5_FOKI[1:], "r1")]
    bottom = [seg("NNDC cohesive end", "NNDC", placeholder=True, length_bp=4),
              seg("P5c paired arm", P5C_NNDC_CORE, "r1")]
    sc = Scene()
    sc.strand("P5-FokI", top, label="P5-FokI")
    sc.anneal("P5c-NNDC-FokI", bottom, to="P5-FokI",
              pair=("P5c paired arm", "P5-FokI paired arm"),
              unpaired=("NNDC cohesive end",), label="P5c-NNDC-FokI",
              mod5="/5Phos/")
    sc.ligation_end("P5c-NNDC-FokI", "NNDC cohesive end", kind="cohesive",
                    polarity="5-prime", owner="adapter",
                    text="4-nt degenerate FokI-compatible 5′ overhang")
    sc.mark("P5-FokI", "unpaired 5′ A", "opposite adapter end")
    sc.labels("P5-FokI")
    sc.labels("P5c-NNDC-FokI")
    return AdapterLigationGeometry.known(
        "cohesive", sc,
        statement="The source prints a phosphorylated four-base NNDC 5′ overhang opposite the FokI-generated substrate overhang; P5H/P5Hc variants add the source-listed H:D pair.")


def _final_core(kind: str) -> list[Segment]:
    if kind == "RNA":
        left = [
            seg("P5", il.P5, "p5"),
            seg("i5 carried index", "J" * 8, "cbc", placeholder=True, feature=I5),
            seg("S5", nx.S5, "s5"),
            seg("left mosaic end", nx.ME, "me"),
            seg("transcript cDNA", "X" * 42, placeholder=True),
            seg("VN-anchor complement", "NB", placeholder=True, length_bp=2),
            seg("oligo-dT16 complement", "A" * 16),
            seg("SbfI site", SBFI),
            seg("round-1 barcode", "A" * 4, "cbc", placeholder=True, feature=BC1),
            seg("RNA BC1 exit arm", revcomp(RNA_DISCRIMINATOR + BC1_COMMON)),
        ]
    elif kind == "DNA":
        left = [
            seg("P5", il.P5, "p5"),
            seg("Read-1 arm remainder", il.TRUSEQ_P5_FULL[len(il.P5):], "r1"),
            seg("histone-mark DNA", "X" * 42, placeholder=True),
            seg("right mosaic end reverse complement", nx.ME_RC, "me"),
            seg("NotI site", NOTI),
            seg("round-1 barcode", "A" * 4, "cbc", placeholder=True, feature=BC1),
            seg("DNA BC1 exit arm", revcomp(DNA_DISCRIMINATOR + BC1_COMMON)),
        ]
    else:
        raise ValueError(f"unknown final PairedTag branch {kind!r}")
    return [
        *left,
        seg("round-2 incoming arm", revcomp(ROUND2.barcode_3_arm)),
        seg("round-2 barcode", "B" * 7, "cbc", placeholder=True, feature=BC2),
        seg("round-2 outgoing arm", revcomp(ROUND2.barcode_5_arm)),
        seg("round-3 incoming arm", revcomp(R3_LIGATING_ARM)),
        seg("round-3 barcode", "C" * 7, "cbc", placeholder=True, feature=BC3),
        seg("round-3 UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
        seg("R04 leading A", "A"),
        seg("Index 1 / Read-2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", "I" * 6, "cbc", placeholder=True, feature=I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ]


DNA_LIBRARY = Construct([replace(x,inferred=True) for x in _final_core("DNA")], name="PairedTag DNA library")
RNA_LIBRARY = Construct([replace(x,inferred=True) for x in _final_core("RNA")], name="PairedTag RNA library")
DNA_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
RNA_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
for library, primers in ((DNA_LIBRARY, DNA_PRIMERS), (RNA_LIBRARY, RNA_PRIMERS)):
    errors = sp.verify(library, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError(f"{library.name}: " + "; ".join(errors))


def read2_barcode_starts(library: Construct) -> tuple[int, int, int]:
    """One-based starts for BC3, BC2 and BC1, derived from the final molecule."""
    hit = sp.locate(library, sp.TRUSEQ["R2"])
    if hit is None or hit.strand != "top":
        raise ValueError(f"{library.name}: Read 2 must extend leftward")
    names = ("round-3 barcode", "round-2 barcode", "round-1 barcode")
    return tuple(hit.start - library.span(name)[1] + 1 for name in names)


PUBLISHED_READ2_BARCODE_WINDOWS = ((10, 13), (47, 50), (84, 87))
for library in (DNA_LIBRARY, RNA_LIBRARY):
    starts = read2_barcode_starts(library)
    if not all(lo <= got <= hi for got, (lo, hi) in
               zip(starts, PUBLISHED_READ2_BARCODE_WINDOWS)):
        raise ValueError(f"{library.name}: barcode starts {starts} miss source windows")

FINAL_LIBRARIES = (
    ("Histone-mark DNA library", DNA_LIBRARY, DNA_PRIMERS,
     "INFERRED complete endpoint — zero-phase representative. Read 1 maps antibody-targeted DNA; Read 2 begins with the UMI and traverses BC3, BC2 and BC1.",
     "The DNA branch uses the P5-FokI adapter and P5 Universal primer. Index 1 is six source-defined bases acquired in a seven-cycle index read. Methods give Read 1 as 50 or 100 cycles; the 100-cycle option is drawn."),
    ("Transcriptome library", RNA_LIBRARY, RNA_PRIMERS,
     "INFERRED complete endpoint — zero-phase oligo(dT)-primed representative. Read 1 maps transcript cDNA from the N5/ME end; Read 2 carries the same UMI and three-part nuclear identity. Random-N6-primed molecules replace the VN/poly(A) tract with the six-base priming tract.",
     "The RNA branch receives N5 by tagmentation. The physical i5 primer window is retained; index2acquisition remains unassigned. Methods give Read 1 as 50 or 100 cycles; the 100-cycle option is drawn."),
)
READ_LENGTHS = {
    "Histone-mark DNA library": {"Read 1": 100, "Index 1 (i7)": 7, "Read 2": 100},
    "Transcriptome library": {"Read 1": 100, "Index 1 (i7)": 7, "Read 2": 100},
}
READ_USES = {
    "Histone-mark DNA library": {"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "extract"},
    "Transcriptome library": {"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "extract"},
}


def oligos():
    return (
        f"<b>pMENTs</b> — 5′-/5Phos/{PMENTS}/ddC/-3′",
        f"<b>DNA_#01_RE</b> — 5′-/5Phos/{BC1_COMMON}{DNA_DISCRIMINATOR}[4-nt BC1]{NOTI}{nx.ME}-3′",
        f"<b>RNA_#01_RE</b> — 5′-/5Phos/{BC1_COMMON}{RNA_DISCRIMINATOR}[4-nt BC1]{SBFI}{RNA_DT}-3′",
        f"<b>RNA_#01_NRE</b> — 5′-/5Phos/{BC1_COMMON}{RNA_DISCRIMINATOR}[4-nt BC1]{SBFI}{RNA_N6}-3′",
        f"<b>R02_#01</b> — 5′-/5Phos/{ROUND2.barcode_5_arm}[7-nt BC2]{ROUND2.barcode_3_arm}-3′",
        f"<b>Linker-R02 / Blocker-R02</b> — 5′-{LINKER_R02}-3′ / 5′-{BLOCKER_R02}-3′",
        f"<b>R04_#01</b> — 5′-{R3_READ2_ARM}[10-nt UMI][7-nt BC3][0–3-nt phase]{R3_LIGATING_ARM}-3′",
        f"<b>Linker-R03 / Quencher-R03</b> — 5′-{LINKER_R03}-3′ / 5′-{QUENCHER_R03}-3′",
        f"<b>Anchor-FokI-GH</b> — 5′-{ANCHOR_HANDLE}{ANCHOR_G}*H-3′ (table prints * before the terminal H)",
        f"<b>P5-FokI / P5c-NNDC-FokI</b> — 5′-{P5_FOKI}-3′ / 5′-/5Phos/{P5C_NNDC}-3′",
        f"<b>P5H-FokI / P5Hc-NNDC-FokI</b> — 5′-{P5H_FOKI}-3′ / 5′-/5Phos/{P5HC_NNDC}-3′",
        f"<b>PA-F / PA-R</b> — 5′-{PA_F}-3′ / 5′-{PA_R}-3′",
        f"<b>N5XX</b> — 5′-{il.P5}[8-nt i5]{nx.S5}-3′",
        f"<b>P7XX</b> — 5′-{il.P7}[6-nt i7]{il.TRUSEQ_READ2[:-1]}-3′",
        f"<b>P5 Universal</b> — 5′-{P5_UNIVERSAL[:-1]}*{P5_UNIVERSAL[-1]}-3′ (terminal phosphorothioate)",
    )


def _initial_rows() -> list[Row]:
    chromatin = Scene.duplex([seg("fixed chromatin", "X" * 40,
                                  placeholder=True)], label="fixed chromatin")
    mrna = Scene()
    mrna.strand("RNA", [seg("nuclear RNA", "X" * 28, placeholder=True),
                         seg("poly(A)", "A" * 16)], label="nuclear RNA")
    return [*chromatin.rows(), *mrna.rows()]


PRIMER_TABLE = "Supplementary Table 1, sheet Primer Sequences"


def _n5xx_state():
    """N5XX as printed: P5, an 8-N i5 window, then the 14-nt S5 footprint."""
    state = Construct([Segment("P5", il.P5), Segment("i5", "N" * 8, feature=I5),
                       Segment("S5", nx.S5)])
    if "".join(x.top for x in state.segments) != N5XX:
        raise ValueError("N5XX state no longer matches the declared primer")
    return state


def _p7xx_state():
    """P7XX as printed: P7, a 6-N i7 window, then TruSeq Read 2 minus its last T."""
    state = Construct([Segment("P7", il.P7), Segment("i7", "N" * 6, feature=I7),
                       Segment("TruSeq Read 2 (minus 3' T)", il.TRUSEQ_READ2[:-1])])
    if "".join(x.top for x in state.segments) != P7XX:
        raise ValueError("P7XX state no longer matches the declared primer")
    return state


def population(rounds=3,*,tailed=False,extended=False,extra_end=""):
    sc=Scene()
    for kind in ("DNA","RNA","RNA-N6"):
        parts=list(_indexed_construct(kind,rounds))
        if tailed or extended:parts.append(role("variable TdT/Anchor native end extent unassigned",inferred=True))
        if extra_end:parts.append(role(extra_end,inferred=True))
        sc.strand(kind,parts,label=kind+" independent identifier-bearing strand class")
        if extended:sc.strand(kind+" newcopies",[role("new complementary/Anchor-extension strands; complete native assembly unassigned",inferred=True)],label="independent copied "+kind+" population")
    return sc

def branch_population(kind,end):
    sc=Scene();sc.strand("target",[*list(_indexed_construct(kind,3)),role(end,inferred=True)],label=kind+" target lineage; native end extent unassigned")
    if kind=="RNA":sc.strand("N6target",[*list(_indexed_construct("RNA-N6",3)),role(end,inferred=True)],label="independent random-hexamer RNA target lineage")
    sc.strand("cleaved background",[role("opposite-modality cleavage fragments; quantitative recovery unassigned")],label="independent carried cleavage population")
    return sc

def cleanup(g,key,template,recipe,library=False):
    rule=g.input(key+"_rule",recipe,MaterialKind.SELECTION)
    out=g.product(key+"_clean","source recovered "+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(key+"_waste","physically excluded bead/column/wash fraction",MaterialKind.WASTE)
    g.react(key,"Recover source-selected fraction",StepKind.CLEANUP,((template,InputRole.TEMPLATE),(rule,InputRole.SELECTION)),(out,waste),changes_molecule=False);g.discard(waste,"Source physical bead/column/size transfer, not cleavage/nonamplification/computational filtering")
    return out

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nuclei=g.input("nuclei","antibody-stained permeabilized nuclei; independent RNA/chromatin",MaterialKind.BIOLOGICAL,state=_initial_rows(),primary=True)
    trans=g.input("bc1_transposome","source BC1pA-Tn5 duplex; physical nontransferredpMENTsddC retained",MaterialKind.ADAPTER,state=tn5_loading_scene())
    tn5=g.input("tn5","source loadedpA-Tn5;37°C60min Mgactivation/EDTAstop",MaterialKind.CATALYST)
    ts=Scene();ts.strand("DNA",list(_indexed_construct("DNA",1)),label="transferredDNA strand; native gap geometry unassigned");ts.strand("RNA",[role("nuclearRNA")],label="independent unconvertedRNA")
    tagged=g.product("tagged","nuclei retaining taggedDNA andRNA",MaterialKind.MOLECULAR,state=ts)
    other=g.product("tn5_other","otherpA-Tn5 endclasses",MaterialKind.MOLECULAR)
    g.react("tagment","Tether and transpose at antibody-bound histone mark",StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),(trans,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Amplification eligibility is not physical removal; recovery through source nuclei washes unassigned."))
    dt=g.input("dt_primer","sourceBC1dT16VN primer",MaterialKind.OLIGO,state=Construct(_bc1_segments("RNA")[:-1]))
    n6=g.input("n6_primer","sourceBC1randomN6 primer",MaterialKind.OLIGO,state=Construct(_bc1_segments("RNA-N6")[:-1]))
    rt=g.input("rt","MaximaHMinusEP0751;source50°C/rampcycling",MaterialKind.CATALYST)
    copied=g.product("copied","DNA plus independently dT/N6-primed single-strandcDNA",MaterialKind.MOLECULAR,state=population(1))
    g.react("reverse_transcribe","Copy nuclearRNA with bothBC1primer families",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(dt,InputRole.STRUCTURAL),(n6,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
    current=cleanup(g,"rt_recovery",copied,"sourcepool/Triton/1000g10min4°C nuclei recovery")
    ligase=g.input("barcode_ligase","T4DNAligaseM0202L;37°C30min",MaterialKind.CATALYST)
    for number in (2,3):
        reagent=Scene()
        reagent.strand("barcode",ROUND2.product_segments() if number==2 else _r3_segments(),label="sourceBC"+str(number),mod5="/5Phos/" if number==2 else "")
        reagent.strand("splint",[Segment("linker",LINKER_R02 if number==2 else LINKER_R03)],label="independent linker")
        bc=g.input("bc"+str(number),"published barcode andlinker reagents",MaterialKind.OLIGO,state=reagent)
        joined=g.product("round"+str(number),"mixedDNA/cDNA afterBC"+str(number),MaterialKind.MOLECULAR,state=population(number))
        g.react("ligate_bc"+str(number),"Ligate sourceBC"+str(number),StepKind.LIGATION,((current,InputRole.TEMPLATE),(bc,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(joined,),geometry=_round2_geometry() if number==2 else _round3_geometry())
        stop=g.input("stop"+str(number),"sourceBlockerR02" if number==2 else "sourceQuencherR03/R04terminator name withEDTA",MaterialKind.OLIGO,state=BLOCKER_R02 if number==2 else QUENCHER_R03)
        stopped=g.product("stopped"+str(number),"indexed mixture after sourceblocker/quencher",MaterialKind.MOLECULAR,state=joined.state)
        g.react("stop_bc"+str(number),"Add sourceblocker/quencher after30min",StepKind.HYBRIDIZATION,((joined,InputRole.TEMPLATE),(stop,InputRole.STRUCTURAL)),(stopped,));current=stopped
    rule=g.input("sublibrary_rule","source count/divide2k–5k nuclei;3.5k representative",MaterialKind.CONDITION)
    chosen=g.product("chosen","chosen mixed nuclei sublibrary",MaterialKind.MOLECULAR,state=current.state)
    rest=g.product("other_sublibraries","other nuclei sublibraries",MaterialKind.MOLECULAR,state=current.state)
    g.react("sublibraries","Divide barcoded nuclei into sublibraries",StepKind.SPLIT,((current,InputRole.TEMPLATE),(rule,InputRole.CONDITION)),(chosen,rest),changes_molecule=False);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,"Unselected sublibrary handling unassigned; not computationally excluded material."))
    lys=g.input("lysis","NaCl/SDS/PK55°C2h",MaterialKind.CONDITION)
    released=g.product("released","released indexed mixedDNA/cDNA",MaterialKind.MOLECULAR,state=chosen.state)
    g.react("lyse","Lyse selected sublibrary",StepKind.DEPROTEINIZATION,((chosen,InputRole.TEMPLATE),(lys,InputRole.CONDITION)),(released,))
    released=cleanup(g,"lysis_cleanup",released,"1×SPRI;12.5µlwater")
    heat=g.input("denature","95°C5min,ice5min",MaterialKind.CONDITION)
    den=g.product("denatured","denatured sourceDNA/cDNA strands",MaterialKind.MOLECULAR,state=released.state)
    g.react("denature_step","Denature beforeTdT",StepKind.THERMAL_TREATMENT,((released,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),(den,))
    tdt=g.input("tdt","TdT+dCTP;37°C30min,75°C20min",MaterialKind.CATALYST)
    tailed=g.product("tailed","variablepolyC-tailed DNA andcDNA",MaterialKind.MOLECULAR,state=population(tailed=True))
    g.react("tail","Add variable3′polyC tails",StepKind.TAILING,((den,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
    anchor=g.input("anchor","sourceAnchor-FokI-GH;terminal*H retained asprinted",MaterialKind.OLIGO,state=ANCHOR_WRITTEN,source=PRIMER_TABLE)
    pol=g.input("extension_polymerase","KAPAHiFiHS;16sourceAnchorcycles",MaterialKind.CATALYST)
    extended=g.product("extended","Anchor-extended mixed targetclasses;nativegeometryunassigned",MaterialKind.MOLECULAR,state=population(extended=True))
    g.react("extend_anchor","Extend from sourceAnchor",StepKind.EXTENSION,((tailed,InputRole.TEMPLATE),(anchor,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(extended,))
    paf=g.input("pa_f","sourcePA-F",MaterialKind.OLIGO,state=PA_F)
    par=g.input("pa_r","sourcePA-R",MaterialKind.OLIGO,state=PA_R)
    preamp=g.product("preamp","pre-amplified mixedDNA/cDNA",MaterialKind.MOLECULAR,state=population(extended=True))
    g.react("preamp_pcr","Preamplify PA-F/PA-R11cycles",StepKind.PCR,((extended,InputRole.TEMPLATE),(paf,InputRole.STRUCTURAL),(par,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(extended.state,reason="ActualAnchor-extension mixture, variabletail/native ends unassigned. No fixed9C or expectedfinalduplex supplies a precursor."))
    preamp=cleanup(g,"preamp_cleanup",preamp,"0.2×+0.65×SPRI double-size;35µlwater")
    dna=g.product("dna_branch","17µlDNA preparation aliquot; mixedmodalitiesretained",MaterialKind.MOLECULAR,state=preamp.state)
    rna=g.product("rna_branch","17µlRNA preparation aliquot; mixedmodalitiesretained",MaterialKind.MOLECULAR,state=preamp.state)
    remain=g.product("split_remainder","remaining1µl preamp product",MaterialKind.MOLECULAR,state=preamp.state)
    g.react("split","Divide mixedpreamp product into sourcealiquots",StepKind.SPLIT,((preamp,InputRole.TEMPLATE),),(dna,rna,remain),changes_molecule=False);g.mark_terminal(Terminal(remain.id,Disposition.UNRESOLVED,"35µl eluate minus two17µl inputs; handlingremaining1µl unassigned."))
    for prefix,template,kind,final,primers,heading in (("dna",dna,"DNA",DNA_LIBRARY,DNA_PRIMERS,"Histone-mark DNA library"),("rna",rna,"RNA",RNA_LIBRARY,RNA_PRIMERS,"Transcriptome library")):
        enzyme=g.input(prefix+"_restriction","SbfI-HF+FokI" if kind=="DNA" else "NotI-HF",MaterialKind.CATALYST)
        cut=g.product(prefix+"_cut","digested mixed"+prefix+" preparation",MaterialKind.MOLECULAR,state=branch_population(kind,"native restriction-produced end/tail extent unassigned"))
        g.react("digest_"+prefix,"Digest source aliquot37°C60min",StepKind.DIGESTION,((template,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(cut,))
        cut=cleanup(g,prefix+"_digest_cleanup",cut,"1.25×SPRI;10µlelution")
        if kind=="DNA":
            adapter=g.input("p5_adapter","sourceP5-FokI/P5c-NNDC duplex; literal variantslisted",MaterialKind.ADAPTER,state=p5_foki_geometry().scenes[0])
            lig=g.input("adapter_ligase","T4DNAligase;4°C10min,10°C15min,16°C15min,25°C45min",MaterialKind.CATALYST)
            tagged=g.product("dna_adapted","P5-adapter-ligated mixedDNApreparation",MaterialKind.MOLECULAR,state=branch_population(kind,"newP5-associated ligated end; nativejunction geometry unassigned"))
            g.react("ligate_p5","Ligate sourceP5adapter",StepKind.LIGATION,((cut,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),(tagged,),geometry=p5_foki_geometry())
            tagged=cleanup(g,"dna_ligation_cleanup",tagged,"1.25×SPRI;30µlwater")
        else:
            adapter=g.input("n5_transposome","sourceTn5-N5 loadedadapter; complete nativecompositionunassigned",MaterialKind.ADAPTER)
            enzyme=g.input("rna_tn5","sourceTn5-N5;37°C30min",MaterialKind.CATALYST)
            tagged=g.product("rna_tagged","pre-indexN5taggedRNA targetclasses plusbackground",MaterialKind.MOLECULAR,state=branch_population(kind,"newN5end; native strand/gaps/fragmentextent unassigned"))
            other=g.product("rna_other_ends","otherRNA-Tn5endclasses",MaterialKind.MOLECULAR)
            g.react("tagment_rna","Tagment mixedRNA-preparation aliquot",StepKind.TAGMENTATION,((cut,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(tagged,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Nonamplification is not physicaldiscard; sourcecleanup recovery unassigned."))
            tagged=cleanup(g,"rna_tagment_cleanup",tagged,"QIAquickPCRpurification;30µl0.1×EB")
        f=g.input(prefix+"_forward","sourceP5Universal (*terminalbond)" if kind=="DNA" else "sourceN5XX8Nindexprimer",MaterialKind.OLIGO,state=Construct([Segment("P5Universal",P5_UNIVERSAL,note="source terminal*mark; exact bond interpretation unassigned")]) if kind=="DNA" else _n5xx_state())
        r=g.input(prefix+"_reverse","sourceP7XX6Nindexprimer;32ntR2arm",MaterialKind.OLIGO,state=_p7xx_state())
        q5=g.input(prefix+"_polymerase","NEBQ5;72°C5min then source8–13additionalStep2repeats to10nM",MaterialKind.CATALYST)
        amp=g.product(prefix+"_amplified","conditional"+prefix+" libraryendpoint component",MaterialKind.MOLECULAR,state=final,inferred=True)
        g.react(prefix+"_pcr","INFERRED native assembly — sourceindexedlibraryPCR",StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(tagged.state,reason="Actualpre-index mixedlibraryprep retained; complete nativegap/ligation/first-cycleassembly is unassigned, not an expectedfinaltemplate."))
        lib=cleanup(g,prefix+"_final_cleanup",amp,"0.9×SPRI",library=True)
        g.sequence(lib,SequencingPlan(library=final,primers=primers,acquired_roles=tuple(p.role for p in primers),read_lengths=READ_LENGTHS[heading],read_uses=READ_USES[heading]))
    g.annotate(SourceDiscrepancy("reverse_transcribe","Oligo(dT) length differs between Methods and Supplementary Table 1",
        (SourceClaim("Online Methods, Reverse transcription (PMC7954905 JATS text line 138)","2.5 µM barcoded T15 primer"),
         SourceClaim("Supplementary Table 1, sheet BC Plate #01, rows B1–B12 (RNA_#01_RE–RNA_#12_RE)","SbfI site followed by 16 T and VN")),
        "The drawn RT primer uses the Supplementary Table 1 sequence (T16VN).",
        "Only the table prints the base-level primer; the Methods value is a name in prose.",
        "Whether 'T15' is a loose name or a different synthesized primer is not stated; the step-by-step Protocol Exchange protocol was not inspected."))
    g.annotate(SourceDiscrepancy("anchor","Anchor oligo name and terminal mark",
        (SourceClaim("Online Methods, Pre-amplification (JATS text line 146)","Anchor-FokI-GSH-Oligo"),
         SourceClaim("Supplementary Table 1, sheet Primer Sequences, row Anchor-FokI-GH","…GAAGGATGTGGGGGGGGG*H")),
        "The table sequence is drawn; '*' is shown as printed before the terminal H.",
        "The table is the only base-level source; the same table writes the P5 Universal primer with a final '*'.",
        "The table does not define '*' or the 'S' in the Methods name; a phosphorothioate reading is an inference, not a source statement."))
    g.annotate(Highlight("reverse_transcribe","RNA-template fate unassigned","RT creates single-strandcDNA; subsequent RNArecovery/degradation is not physicaldiscard in thismodel."))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    return g.require_publishable()

SHOW_READ_CONTENT=False
READ_SPECS={name:(sp.read_spec("Read 1",100,"targetDNA or transcriptcDNA;50cycle publishedalternative retained innarrative","include"),sp.read_spec("Index 1 (i7)",7,"source6ntindex plus adjacentadapterbase (derived seventhcycle)","demultiplex"),sp.read_spec("Read 2",100,"N10UMI andBC3/BC2/BC1 chain","extract"),sp.read_spec("Index 2 (i5)","unassigned","RNAphysical8ntprimerwindow; acquisitionunassigned (not applicable toDNAprimerpair)","undetermined")) for name,*_ in FINAL_LIBRARIES}
READ_SPECS_NOTE="Published100+7+100 example retained;50+7+100 is a separate positive source option. Index2acquisition remainsunassigned."
SOURCE_OLIGO_EXCEPTIONS={"bc1_transposome":"Physical pMENTs terminalddC is a nucleotide; scraper._letters strips it withoutC. Direct wholemodifiedpMENTs and transferredstrand source regressions cover thisinput; no alteredsequence supplied."}
REVIEWED_TARGETS=("nuclei","bc1_transposome","tn5","tagment","dt_primer","n6_primer","rt","reverse_transcribe","rt_recovery_rule","rt_recovery","barcode_ligase","bc2","ligate_bc2","stop2","stop_bc2","bc3","ligate_bc3","stop3","stop_bc3","sublibrary_rule","sublibraries","lysis","lyse","lysis_cleanup_rule","lysis_cleanup","denature","denature_step","tdt","tail","anchor","extension_polymerase","extend_anchor","pa_f","pa_r","preamp_pcr","preamp_cleanup_rule","preamp_cleanup","split","dna_restriction","digest_dna","dna_digest_cleanup_rule","dna_digest_cleanup","rna_restriction","digest_rna","rna_digest_cleanup_rule","rna_digest_cleanup","p5_adapter","adapter_ligase","ligate_p5","dna_ligation_cleanup_rule","dna_ligation_cleanup","n5_transposome","rna_tn5","tagment_rna","rna_tagment_cleanup_rule","rna_tagment_cleanup","dna_forward","dna_reverse","dna_polymerase","dna_pcr","dna_final_cleanup_rule","dna_final_cleanup","dna_final_cleanup_clean","rna_forward","rna_reverse","rna_polymerase","rna_pcr","rna_final_cleanup_rule","rna_final_cleanup","rna_final_cleanup_clean")
