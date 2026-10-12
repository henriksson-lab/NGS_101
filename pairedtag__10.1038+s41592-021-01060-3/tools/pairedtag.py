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
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      duplex_rows, feature, revcomp)
from multimodal_spatial import HybridizationLigationRound


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

ROUND2 = HybridizationLigationRound(
    "round 2", BC1_COMMON, "GTGCGAACTCAGACC", "ATCCACGTGCTTGAG",
    LINKER_R02,
    feature("pairedtag_cell_bc2", "cell_barcode", "whitelist",
            group="pairedtag_cell_barcode", part="round 2",
            whitelist="Zhu et al. Supplementary Table 1, BC Plate #02"),
    phosphorylated_5=True, barcode_length=7)

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
    bottom = [seg("pMENTs mosaic-end complement", PMENTS, "me")]
    sc = Scene()
    sc.strand("barcoded top", top, label="DNA_#01_RE", mod5="/5Phos/")
    sc.anneal("pMENTs", bottom, to="barcoded top",
              pair=("pMENTs mosaic-end complement", "mosaic end"),
              label="pMENTs", mod5="/5Phos/", mod3="ddC")
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


def _tailed_scene(kind: str) -> Scene:
    parts = [*list(_indexed_construct(kind, 3)), seg("TdT poly(C) tail", "C" * 9)]
    sc = Scene()
    sc.strand(kind.lower(), parts, label=f"{kind} after TdT")
    sc.mark(kind.lower(), "TdT poly(C) tail", "3′ poly(C) tail")
    sc.labels(kind.lower())
    return sc


def _anchor_scene(kind: str) -> Scene:
    parts = [*list(_indexed_construct(kind, 3)), seg("TdT poly(C) tail", "C" * 9)]
    anchor = [seg("anchor / PA-R handle", ANCHOR_HANDLE, "tso"),
              seg("anchor G9", ANCHOR_G),
              seg("anchor H", "H", placeholder=True)]
    sc = Scene()
    sc.strand(kind.lower(), parts, label=f"{kind} poly(C)-tailed strand")
    sc.anneal("Anchor-FokI-GH", anchor, to=kind.lower(),
              pair=("anchor G9", "TdT poly(C) tail"),
              unpaired=("anchor / PA-R handle", "anchor H"),
              label="Anchor-FokI-GH")
    sc.arrow("Anchor-FokI-GH", "linear extension from the anchored 3′ end")
    sc.labels("Anchor-FokI-GH")
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


DNA_LIBRARY = Construct(_final_core("DNA"), name="PairedTag DNA library")
RNA_LIBRARY = Construct(_final_core("RNA"), name="PairedTag RNA library")
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
     "Zero-phase representative. Read 1 maps antibody-targeted DNA; Read 2 begins with the UMI and traverses BC3, BC2 and BC1.",
     "The DNA branch uses the P5-FokI adapter and P5 Universal primer. Index 1 is six source-defined bases acquired in a seven-cycle index read."),
    ("Transcriptome library", RNA_LIBRARY, RNA_PRIMERS,
     "Zero-phase oligo(dT)-primed representative. Read 1 maps transcript cDNA from the N5/ME end; Read 2 carries the same UMI and three-part nuclear identity. Random-N6-primed molecules replace the VN/poly(A) tract with the six-base priming tract.",
     "The RNA branch receives N5 by tagmentation. The carried i5 is not acquired in the published single-index run."),
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
        f"<b>Anchor-FokI-GH</b> — 5′-{ANCHOR_WRITTEN}-3′ (penultimate bond phosphorothioate)",
        f"<b>P5-FokI / P5c-NNDC-FokI</b> — 5′-{P5_FOKI}-3′ / 5′-/5Phos/{P5C_NNDC}-3′",
        f"<b>P5H-FokI / P5Hc-NNDC-FokI</b> — 5′-{P5H_FOKI}-3′ / 5′-/5Phos/{P5HC_NNDC}-3′",
        f"<b>PA-F / PA-R</b> — 5′-{PA_F}-3′ / 5′-{PA_R}-3′",
        f"<b>N5XX</b> — 5′-{il.P5}[8-nt i5]{nx.S5}-3′",
        f"<b>P7XX</b> — 5′-{il.P7}[6-nt i7]{il.TRUSEQ_READ2[:-1]}-3′",
        f"<b>P5 Universal</b> — 5′-{P5_UNIVERSAL}-3′",
    )


def _initial_rows() -> list[Row]:
    chromatin = Scene.duplex([seg("fixed chromatin", "X" * 40,
                                  placeholder=True)], label="fixed chromatin")
    mrna = Scene()
    mrna.strand("RNA", [seg("nuclear RNA", "X" * 28, placeholder=True),
                         seg("poly(A)", "A" * 16)], label="nuclear RNA")
    return [*chromatin.rows(), *mrna.rows()]


def workflow() -> Workflow:
    wf = Workflow(MolecularState("Antibody-stained permeabilized nuclei",
                                 tuple(_initial_rows())))
    wf.react("Anneal the BC1 transposon and load pA–Tn5",
             tn5_loading_scene().rows(), name="BC1-loaded pA–Tn5",
             note="DNA_#01_RE carries the exact common arm, four-base BC1, NotI site and mosaic end; pMENTs is 5′-phosphorylated and 3′-ddC blocked.")
    wf.react("Tether and transpose at the antibody-bound histone mark",
             _indexed_scene("DNA", 1).rows(), name="BC1-tagged histone-mark DNA")
    wf.react("Reverse-transcribe RNA with both BC1 primer families",
             [*rt_scene(False).rows(), *rt_scene(True).rows(),
              *_both_modalities(1)], name="BC1-tagged RNA and DNA",
             note="Matched four-base BC1 sequences occur in the oligo(dT)16-VN, random-N6 and transposome oligos.")
    wf.ligate_adapter("Splint-ligate BC2", _round2_geometry(),
                      _both_modalities(2), name="After barcode round 2",
                      note="Blocker-R02 sequesters the linker after the first ligation.")
    wf.ligate_adapter("Splint-ligate BC3, its UMI and phase block",
                      _round3_geometry(), _both_modalities(3),
                      name="Three-round indexed molecules",
                      note="The zero-phase A1 member is drawn; other wells add one to three phase bases before the common ligating arm.")
    wf.react("Lyse nuclei and add a 3′ poly(C) tail with TdT",
             [*_tailed_scene("RNA").rows(), *_tailed_scene("DNA").rows()],
             name="Poly(C)-tailed indexed molecules")
    wf.react("Anneal Anchor-FokI-GH and linearly extend",
             [*_anchor_scene("RNA").rows(), *_anchor_scene("DNA").rows()],
             name="Anchor-extended products",
             note="The source oligo contains the PA-R handle, FokI element, G9 tract and terminal H.")
    wf.react("Preamplify with PA-F/PA-R and split the two modalities",
             [Row(chunks=[("DNA aliquot: SbfI + FokI expose the DNA finishing end", None, False)]),
              Row(chunks=[("RNA aliquot: NotI exposes the RNA product for N5-Tn5", None, False)])],
             name="Restriction-separated DNA and RNA branches")
    wf.ligate_adapter("Ligate the P5-FokI duplex to the DNA branch",
                      p5_foki_geometry(),
                      [*duplex_rows(DNA_LIBRARY, label="DNA library"),
                       Row(chunks=[("RNA branch awaits N5 tagmentation", None, False)])],
                      name="P5-ligated DNA branch and digested RNA branch",
                      note="The adapter mixture combines P5/P5c and P5H/P5Hc at 1:3; the unphased representative shows the common NNDC cohesive end.")
    wf.react("Tagment the RNA branch with Tn5-N5 and index both libraries",
             [*duplex_rows(DNA_LIBRARY, label="DNA library"),
              *duplex_rows(RNA_LIBRARY, label="RNA library")],
             name="PairedTag sequencing libraries",
             note="The source run is PE 50 or 100 + 7 + 100 (Read 1 + Index 1 + Read 2).")
    return wf
