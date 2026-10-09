"""Molecular construct model for PETRI-seq (Blattman et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature

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
