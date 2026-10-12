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
from multimodal_spatial import HybridizationLigationRound

TITLE = "Paired-seq — combinatorial RNA and chromatin accessibility"
NOTES = "01_paired-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41594-019-0323-x">'
          'Zhu et al., <i>Nature Structural &amp; Molecular Biology</i> (2019)</a>.')
SUMMARY = ("Permeabilized nuclei receive a three-base first barcode through indexed "
           "RNA reverse transcription or chromatin tagmentation. Three splinted "
           "ligation rounds add BC2–BC4 and a UMI; TdT/anchor amplification is then "
           "split by restriction digest into RNA and ATAC libraries.")
CAVEAT = ("Round-4 oligos contain a 0–3-nt phase block. Drawings use a published "
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

ROUND2 = HybridizationLigationRound(
    "round 2", BC1_COMMON_ARM, "GTGCGAACTCAGACC", "ATCCACGTGCTTGAG",
    "CGAATGCTCTGGCCTCTCAAGCACGTGGAT", BC2, phosphorylated_5=True,
    barcode_length=7)
ROUND3 = HybridizationLigationRound(
    "round 3", ROUND2.barcode_5_arm, "CATCGGCGTACGACT", "GGATTCGAGGAGCGT",
    "GGTCTGAGTTCGCACACGCTCCTCGAATCC", BC3, phosphorylated_5=True,
    barcode_length=7)

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


def _tailed_scene(kind: str) -> Scene:
    parts = [*list(_indexed_construct(kind, 4)),
             _seg("TdT poly(C) tail", "C" * 9)]
    sc = Scene()
    sc.strand(kind.lower(), parts, label=f"{kind} after TdT")
    sc.mark(kind.lower(), "TdT poly(C) tail", "3′ poly(C) tail")
    sc.labels(kind.lower())
    return sc


def _anchor_scene(kind: str) -> Scene:
    parts = [*list(_indexed_construct(kind, 4)),
             _seg("TdT poly(C) tail", "C" * 9)]
    anchor = [_seg("anchor handle", ANCHOR_HANDLE, "tso"),
              _seg("anchor G tract", ANCHOR_G),
              _seg("anchoring H", "H", placeholder=True)]
    sc = Scene()
    sc.strand(kind.lower(), parts, label=f"{kind} poly(C)-tailed strand")
    sc.anneal("anchor", anchor, to=kind.lower(),
              pair=("anchor G tract", "TdT poly(C) tail"),
              label="Anchor oligo", unpaired=("anchor handle", "anchoring H"))
    sc.arrow("anchor", "linear extension from the anchored 3′ end")
    sc.labels("anchor")
    return sc


def _final_core(kind: str) -> list[Segment]:
    if kind == "RNA":
        branch = [_seg("SbfI site", SBFI),
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
        _seg("dA junction", "A"),
        _seg("TruSeq Read-2 arm", il.CANONICAL_ARM_READ2, "r2"),
        _seg("i7 reverse complement", "I" * 6, "cbc", placeholder=True,
             molecular_feature=I7_FEATURE),
        _seg("P7 reverse complement", il.P7_RC, "p7"),
    ]


PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def _library(kind: str) -> Construct:
    lib = Construct(_final_core(kind), name=f"Paired-seq {kind} library")
    problems = sp.verify(lib, PRIMERS,
                         required_roles=tuple(p.role for p in PRIMERS))
    if problems:
        raise ValueError(f"invalid Paired-seq {kind} library: " + "; ".join(problems))
    return lib


RNA = _library("RNA")
ATAC = _library("ATAC")
FINAL_LIBRARIES = (
    ("RNA library", RNA, PRIMERS,
     "Zero-phase representative. Read 1 enters cDNA; Read 2 begins with UMI and the four-part cell barcode chain.",
     "The published mixed primer set is derived from the completed library."),
    ("Chromatin-accessibility library", ATAC, PRIMERS,
     "Zero-phase representative. Read 1 enters accessible DNA; Read 2 begins with UMI and the matched four-part cell barcode chain.",
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
        f"<b>P7XX</b> — 5′-{il.P7}[6-nt i7]{il.TRUSEQ_READ2}-3′",
    )


def _initial_rows() -> list[Row]:
    chromatin = Scene.duplex([_seg("fixed chromatin", "X" * 40,
                                  placeholder=True)], label="fixed chromatin")
    mrna = Scene()
    mrna.strand("mRNA", [_seg("transcript", "X" * 28, placeholder=True),
                         _seg("poly(A)", "A" * 16)], label="poly(A) RNA")
    return [*chromatin.rows(), *mrna.rows()]


def workflow() -> Workflow:
    wf = Workflow(MolecularState("Permeabilized nuclei", tuple(_initial_rows())))
    wf.react(
        "Transpose chromatin with BC1-indexed Tn5",
        _indexed_scene("ATAC", 1).rows(), name="BC1-tagged accessible DNA",
        note="The transferred strand contains the exact common arm, 3-nt BC1, NotI site and mosaic end; pMENTs is 5′-phosphorylated and 3′-ddC blocked.")
    wf.react(
        "Reverse-transcribe RNA with the BC1-indexed oligo-dT/VN primer",
        _both_modalities(1), name="BC1-tagged RNA and ATAC molecules",
        note="The RNA primer uses the same 15-nt common arm and logical BC1, but carries SbfI plus 16T-VN instead of NotI plus mosaic end.")
    for number in (2, 3):
        wf.ligate_adapter(
            f"Splint-ligate barcode round {number}", _geometry(number),
            _both_modalities(number), name=f"After barcode round {number}",
            note=f"Linker-R0{number} aligns the published barcode oligo with the preceding 5′ arm; ** marks the sealed nick.")
    wf.ligate_adapter(
        "Splint-ligate BC4, its UMI and phase block", _round4_geometry(),
        _both_modalities(4), name="Four-round indexed molecules",
        note="The displayed zero-phase member carries the exact Read-2 arm, 10-nt UMI and 7-nt BC4 before the final ligating arm.")
    wf.react(
        "Lyse nuclei and add a 3′ poly(C) tail with TdT",
        [*_tailed_scene("RNA").rows(), *_tailed_scene("ATAC").rows()],
        name="Poly(C)-tailed indexed molecules",
        note="TdT creates the substrate for the G-rich Anchor oligo on both analytes.")
    wf.react(
        "Anneal Anchor and perform linear extension",
        [*_anchor_scene("RNA").rows(), *_anchor_scene("ATAC").rows()],
        name="Anchor-extended products",
        note="The nine-G tract is paired to a representative nine-C tail; terminal H is the degenerate 3′ anchor base.")
    wf.react(
        "Preamplify with PA-F and PA-R, then split by SbfI or NotI",
        [Row(chunks=[("RNA aliquot: NotI removes the DNA architecture", "r2", False)]),
         Row(chunks=[("ATAC aliquot: SbfI removes the RNA architecture", "me", False)])],
        name="Modality-separated amplified products",
        note="The mutually exclusive restriction sites introduced with BC1 enable physical branch selection.")
    wf.react(
        "Retagment with Tn5-P5 and PCR-complete N5XX/P7XX",
        [*duplex_rows(RNA, label="RNA library"),
         *duplex_rows(ATAC, label="ATAC library")],
        name="Paired-seq sequencing libraries",
        note="The result combines a Nextera Read-1 end with a TruSeq Read-2 end. The defining run acquired 53 + 7 + 130 cycles (R1 + i7 + R2).")
    return wf
