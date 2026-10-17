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
from multimodal_spatial import HybridizationLigationRound, finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

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
N5XX = il.P5 + "N" * 8 + nx.S5
P7XX = il.P7 + "N" * 6 + il.TRUSEQ_READ2
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


def _single_strand(name: str, bases: str) -> Scene:
    sc = Scene()
    sc.strand(name, [_seg(name, bases)], label=name)
    return sc


def _preamp_template(kind: str) -> Construct:
    """Exact anchor-extended substrate used by the shared PA-F/PA-R PCR."""
    return Construct([*list(_indexed_construct(kind, 4)),
                      _seg("TdT poly(C) tail", "C" * 9),
                      _seg("anchor handle complement", revcomp(PA_R))],
                     name=f"Paired-seq {kind} anchor-extended preamplification template")


def protocol_graph():
    """Directed provenance for the paired RNA/ATAC chemistry.

    The graph keeps nuclei intact through four barcode coordinates, then makes the
    restriction-selected RNA and ATAC branches explicit.  Every tagmentation exposes
    its non-library products and both final libraries carry independently checked
    sequencing plans.
    """
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = g.input(
        "nuclei", "permeabilized nuclei containing chromatin and poly(A) RNA",
        MaterialKind.BIOLOGICAL, state=_initial_rows(), primary=True)

    transposome = g.input(
        "bc1_transposome", "BC1-indexed Tn5 (DNA_#01-08_RE annealed to pMENTs), 37 °C 30 min, EDTA stop",
        MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "unloaded Tn5 transposase (loaded 30 min RT)", MaterialKind.CATALYST)
    tagged = g.product(
        "tagged", "nuclei retaining BC1-tagged chromatin and RNA",
        MaterialKind.MOLECULAR,
        state=[*_indexed_scene("ATAC", 1).rows(),
               Row(chunks=[("poly(A) RNA remains available for reverse transcription",
                            None, False)])])
    tn5_waste = g.product(
        "tn5_waste", "non-library chromatin tagmentation products",
        MaterialKind.WASTE)
    g.react(
        "tagment_chromatin", "Transpose chromatin with BC1-indexed Tn5",
        StepKind.TAGMENTATION,
        ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
         (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    g.discard(tn5_waste, "does not retain an amplifiable BC1-tagged chromatin fragment")

    rt_primer = g.input(
        "rt_primer", "matched RNA_#01-08_RE oligo-dT/VN primer (same well order as Tn5 BC1)",
        MaterialKind.OLIGO)
    rt = g.input("rt", "Maxima reverse transcriptase, multi-temperature cycling program", MaterialKind.CATALYST)
    copied = g.product(
        "copied", "nuclei retaining matched BC1-tagged RNA and ATAC molecules",
        MaterialKind.MOLECULAR, state=_both_modalities(1))
    g.react(
        "reverse_transcribe", "Reverse-transcribe RNA with the BC1-indexed primer",
        StepKind.REVERSE_TRANSCRIPTION,
        ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (rt, InputRole.CATALYST)), (copied,))

    ligase = g.input("barcode_ligase", "T4 DNA ligase (NEB M0202L), 37 °C 30 min per round",
                     MaterialKind.CATALYST)
    current = copied
    reviewed = [
        "nuclei", "bc1_transposome", "tn5", "tagment_chromatin",
        "rt_primer", "rt", "reverse_transcribe", "barcode_ligase",
    ]
    for number, geometry in ((2, _geometry(2)), (3, _geometry(3)),
                             (4, _round4_geometry())):
        barcode = g.input(
            f"barcode_{number}", f"round-{number} barcode and splint assembly",
            MaterialKind.OLIGO, state=geometry.scenes[0])
        joined = g.product(
            f"round_{number}", f"RNA and ATAC molecules after barcode round {number}",
            MaterialKind.MOLECULAR, state=_both_modalities(number))
        if number < 4:
            stop = g.input(
                f"blocker_{number}",
                f"Blocker-R0{number} (reverse complement of Linker-R0{number}), added after 30 min",
                MaterialKind.OLIGO,
                state=_single_strand(f"Blocker-R0{number}",
                                     {2: ROUND2, 3: ROUND3}[number].blocking_oligo))
        else:
            stop = g.input(
                "terminator", "R04 Terminator (reverse complement of Linker-R04) with EDTA",
                MaterialKind.OLIGO,
                state=_single_strand("Terminator", revcomp(LINKER_R04)))
        g.react(
            f"ligate_{number}", f"Splint-ligate barcode round {number}",
            StepKind.LIGATION,
            ((current, InputRole.TEMPLATE), (barcode, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST), (stop, InputRole.STRUCTURAL)), (joined,),
            geometry=geometry)
        reviewed.extend((barcode.id, stop.id, f"ligate_{number}"))
        current = joined

    lysis = g.input("lysis_condition", "NaCl/SDS/proteinase K, 55 °C 2 h, then 1× SPRI",
                    MaterialKind.CONDITION)
    released = g.product(
        "released", "released four-round indexed RNA and ATAC molecules",
        MaterialKind.MOLECULAR, state=current.state)
    g.react("lyse", "Lyse nuclei and release indexed molecules",
            StepKind.EXTRACTION,
            ((current, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)),
            (released,))

    tdt = g.input("tdt", "TdT + dCTP after 95 °C denaturation, 37 °C 30 min",
                  MaterialKind.CATALYST)
    tailed = g.product(
        "tailed", "poly(C)-tailed indexed RNA and ATAC molecules",
        MaterialKind.MOLECULAR,
        state=[*_tailed_scene("RNA").rows(), *_tailed_scene("ATAC").rows()])
    g.react("tail", "Add a 3-prime poly(C) tail with TdT", StepKind.TAILING,
            ((released, InputRole.TEMPLATE), (tdt, InputRole.CATALYST)),
            (tailed,))

    anchor = g.input("anchor", "G-rich Anchor oligo", MaterialKind.OLIGO)
    annealed = g.product(
        "anchor_annealed", "Anchor oligo annealed to both poly(C) tails",
        MaterialKind.MOLECULAR,
        state=[*_anchor_scene("RNA").rows(), *_anchor_scene("ATAC").rows()])
    g.react("anneal_anchor", "Anneal the G-rich Anchor oligo",
            StepKind.HYBRIDIZATION,
            ((tailed, InputRole.TEMPLATE), (anchor, InputRole.STRUCTURAL)),
            (annealed,))
    extension_polymerase = g.input(
        "extension_polymerase", "KAPA HiFi HS, 15 linear Anchor-extension cycles",
        MaterialKind.CATALYST)
    extended = g.product(
        "extended", "anchor-extended RNA and ATAC products",
        MaterialKind.MOLECULAR, state=annealed.state)
    g.react("extend_anchor", "Linearly extend from the anchored 3-prime end",
            StepKind.EXTENSION,
            ((annealed, InputRole.TEMPLATE),
             (extension_polymerase, InputRole.CATALYST)), (extended,))

    pa_f = g.input("pa_f", "PA-F preamplification primer", MaterialKind.OLIGO)
    pa_r = g.input("pa_r", "PA-R preamplification primer", MaterialKind.OLIGO)
    preamp_polymerase = g.input(
        "preamp_polymerase", "KAPA HiFi HS, 10 PA-F/PA-R cycles, SPRI double size selection",
        MaterialKind.CATALYST)
    preamp = g.product(
        "preamp", "preamplified Paired-seq RNA and ATAC products",
        MaterialKind.MOLECULAR, state=_both_modalities(4))
    preamp_geometry = PCRGeometry.exact(
        _preamp_template("ATAC"), PCRPrimer("PA-F", PA_F),
        PCRPrimer("PA-R", PA_R),
        statement="The published PA-F and PA-R sites face inward on the anchor-extended product.")
    g.react(
        "preamp_pcr", "Preamplify with PA-F and PA-R", StepKind.PCR,
        ((extended, InputRole.TEMPLATE), (pa_f, InputRole.STRUCTURAL),
         (pa_r, InputRole.STRUCTURAL),
         (preamp_polymerase, InputRole.CATALYST)), (preamp,),
        geometry=preamp_geometry)

    rna = g.product("rna_branch", "RNA finishing aliquot",
                    MaterialKind.MOLECULAR, state=preamp.state)
    atac = g.product("atac_branch", "ATAC finishing aliquot",
                     MaterialKind.MOLECULAR, state=preamp.state)
    g.react("split", "Split preamplified material into RNA and ATAC aliquots",
            StepKind.SPLIT, ((preamp, InputRole.TEMPLATE),), (rna, atac),
            changes_molecule=False)

    rna_re = g.input("rna_restriction", "NotI-HF, 37 °C 60 min, QIAquick", MaterialKind.CATALYST)
    rna_cut = g.product("rna_cut", "NotI-selected RNA branch",
                        MaterialKind.MOLECULAR, state=rna.state)
    g.react("digest_rna", "Digest the RNA aliquot with NotI",
            StepKind.DIGESTION,
            ((rna, InputRole.TEMPLATE), (rna_re, InputRole.CATALYST)),
            (rna_cut,))
    atac_re = g.input("atac_restriction", "SbfI-HF, 37 °C 60 min, QIAquick", MaterialKind.CATALYST)
    atac_cut = g.product("atac_cut", "SbfI-selected ATAC branch",
                         MaterialKind.MOLECULAR, state=atac.state)
    g.react("digest_atac", "Digest the ATAC aliquot with SbfI",
            StepKind.DIGESTION,
            ((atac, InputRole.TEMPLATE), (atac_re, InputRole.CATALYST)),
            (atac_cut,))

    final_transposome = g.input(
        "p5_transposome", "Tn5-P5 finishing transposome (adapter loading not described)",
        MaterialKind.ADAPTER)
    final_tn5 = g.input("finishing_tn5", "finishing Tn5 transposase",
                        MaterialKind.CATALYST)
    branch_templates = []
    for prefix, template, final in (("rna", rna_cut, RNA),
                                    ("atac", atac_cut, ATAC)):
        tagged_branch = g.product(
            prefix + "_tagged", prefix.upper() + " branch after P5 tagmentation",
            MaterialKind.MOLECULAR, state=final)
        waste = g.product(prefix + "_tn5_waste",
                          "non-library " + prefix.upper() + " tagmentation products",
                          MaterialKind.WASTE)
        g.react(
            "tagment_" + prefix,
            "Retagment the " + prefix.upper() + " branch with P5-loaded Tn5",
            StepKind.TAGMENTATION,
            ((template, InputRole.TEMPLATE),
             (final_transposome, InputRole.STRUCTURAL),
             (final_tn5, InputRole.CATALYST)), (tagged_branch, waste))
        g.discard(waste, "does not retain the selected P5-ended insert")
        reviewed.append("tagment_" + prefix)
        branch_templates.append(tagged_branch)

    reviewed.extend((
        "lysis_condition", "lyse", "tdt", "tail", "anchor", "anneal_anchor",
        "extension_polymerase", "extend_anchor", "pa_f", "pa_r",
        "preamp_polymerase", "preamp_pcr", "split", "rna_restriction",
        "digest_rna", "atac_restriction", "digest_atac", "p5_transposome",
        "finishing_tn5"))
    for prefix, template, final, primers, title in (
            ("rna", branch_templates[0], RNA, PRIMERS, "Paired-seq RNA"),
            ("atac", branch_templates[1], ATAC, PRIMERS, "Paired-seq ATAC")):
        geometry = PCRGeometry.exact(
            final, PCRPrimer("N5XX", N5XX, anneal_nt=len(nx.S5)),
            PCRPrimer("P7XX", P7XX),
            statement="The published N5XX and P7XX primer sites face inward on the completed branch.")
        finish_pcr_branch(
            g, prefix + "_final", template, final, primers, reviewed,
            label=title, read_lengths=READ_LENGTHS[
                "RNA library" if prefix == "rna" else "Chromatin-accessibility library"],
            read_uses=READ_USES[
                "RNA library" if prefix == "rna" else "Chromatin-accessibility library"],
            geometry=geometry)

    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
