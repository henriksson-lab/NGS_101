"""Molecular construct model for PETRI-seq (Blattman et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp

UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "fixed")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    if name == "UMI": kw.setdefault("feature", UMI_FEATURE)
    elif name == "i5 index": kw.setdefault("feature", I5_FEATURE)
    elif name == "i7 index read": kw.setdefault("feature", I7_FEATURE)
    return Segment(name=name, top=top, tag=tag, **kw)


R1_PREFIX = "GCCAGA"
R2_PREFIX = "GCTTCGC"
R2_SUFFIX = "CCTCCTAC"
R3_PREFIX = "AGAATACACGACGCTCTTCCGATCT"
R3_SUFFIX = "GGTCCTTG"
LINKER_R2 = "STCTGGCGTAGGAGGW"
LINKER_R3 = "GCGAAGCCAAGGACCW"
I5_OLIGO_INDEX = il.NEBNEXT_I5_SET1["i501"]
I7_OLIGO_INDEX = "TAAGGCGA"
I7_INDEX_READ = "TCGCCTTA"


def barcode(name: str) -> Segment:
    round_name = name.replace("barcode ", "round ")
    return seg(name, "B" * 7, "cbc", placeholder=True,
               feature=feature(name.replace(" ", "_"), "cell_barcode", "combinatorial",
                               group="cell_barcode", part=round_name))


def round1_primer() -> list[Segment]:
    return [seg("round-1 constant", R1_PREFIX), barcode("barcode 1"),
            seg("random hexamer", "N" * 6, placeholder=True)]


def round2_oligo() -> list[Segment]:
    return [seg("round-2 constant", R2_PREFIX), barcode("barcode 2"),
            seg("round-2 junction", R2_SUFFIX)]


def round3_oligo() -> list[Segment]:
    return [seg("round-3 library handle", R3_PREFIX, "r1"),
            seg("UMI", "N" * 7, "umi", placeholder=True), barcode("barcode 3"),
            seg("round-3 junction", R3_SUFFIX)]


def after_round1(insert_nt: int = 28) -> Construct:
    return Construct([*round1_primer(), seg("cDNA", "X" * insert_nt, placeholder=True)],
                     name="PETRI-seq after RT barcode")


def after_round2(insert_nt: int = 28) -> Construct:
    return Construct([*round2_oligo(), *after_round1(insert_nt)],
                     name="PETRI-seq after second barcode")


def after_round3(insert_nt: int = 28) -> Construct:
    return Construct([*round3_oligo(), *after_round2(insert_nt)],
                     name="PETRI-seq after third barcode")


def splint_geometry(round_no: int) -> AdapterLigationGeometry:
    linker = LINKER_R2 if round_no == 2 else LINKER_R3
    rc = revcomp(linker)
    sc = Scene(); sc.strand("splint", [seg(f"round-{round_no} new half", linker[8:], placeholder=True),
                                       seg(f"round-{round_no} old half", linker[:8], placeholder=True)],
                                  label=f"round-{round_no} linker")
    sc.anneal("ligating oligos", [seg(f"round-{round_no} 5-prime arm", rc[:8], placeholder=True),
                                   seg("previous-round 3-prime arm", rc[8:], placeholder=True)],
              to="splint", pair=(f"round-{round_no} 5-prime arm",
                                  f"round-{round_no} new half"),
              label="adjacent barcode oligos", mod5="p",
              unpaired=(f"round-{round_no} 5-prime arm", "previous-round 3-prime arm",
                        f"round-{round_no} new half", f"round-{round_no} old half"))
    sc.ligation_end("ligating oligos", f"round-{round_no} 5-prime arm",
                    kind="splinted", polarity="5-prime", owner="adapter",
                    text="5′-phosphorylated nick partner")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement=f"The PETRI-seq supplement prints the 16-base round-{round_no} linker; its two eight-base halves align the barcode oligos across the nick.")


def workflow() -> Workflow:
    r1, r2, r3 = after_round1(), after_round2(), after_round3()
    wf = Workflow(MolecularState("Round-1 barcoded cDNA", tuple(Scene.duplex(list(r1), label="after RT").rows())))
    wf.ligate_adapter("Ligate barcode round 2", splint_geometry(2),
                      Scene.duplex(list(r2), label="after round 2").rows(),
                      note="The splint spans the round-2/round-1 nick.")
    wf.ligate_adapter("Ligate UMI and barcode round 3", splint_geometry(3),
                      Scene.duplex(list(r3), label="after round 3").rows(),
                      note="The third oligo supplies the UMI, barcode 3 and future Read-1 handle.")
    wf.react("Second-strand synthesis, tagmentation and selective PCR",
             final_scene().rows(), note="NEB i50x and Nextera N7xx select the two required ends.")
    return wf


SEQ_PRIMERS = tuple([sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
                     sp.TRUSEQ["I2"], sp.NEXTERA["R2"]])


def final_library(insert_nt: int = 28) -> Construct:
    """Selected round-3-handle/s7 product after mixed NEB-i5/Nextera-i7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 index", I5_OLIGO_INDEX, "cbc"),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("UMI", "N" * 7, "umi", placeholder=True), barcode("barcode 3"),
        seg("round-3 junction", R3_SUFFIX), seg("round-2 constant", R2_PREFIX),
        barcode("barcode 2"), seg("round-2 junction", R2_SUFFIX),
        seg("round-1 constant", R1_PREFIX), barcode("barcode 1"),
        seg("random-hexamer-derived bases", "N" * 6, placeholder=True),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="PETRI-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final PETRI-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":58})
    by_start={x.cycle_start:(x.cycles,x.feature.label+(
              f" ({x.feature.part})" if x.feature.part else "")) for x in spans}
    return [by_start[1],by_start[8], ("15–29", "round-3/round-2 junction"),
            by_start[30], ("37–50", "round-2/round-1 junction"), by_start[51],
            ("58", "first random-hexamer-derived base")]
