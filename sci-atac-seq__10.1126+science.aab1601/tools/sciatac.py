"""Molecular constructors for the individually published sci-ATAC protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp

SCI15 = "sci-ATAC-seq"
SCI18 = "sci-ATAC-seq (2018)"
SCI3 = "sci-ATAC-seq3"

A_LINK = "TCCACGC"
B_LINK = "CTGTCCCTGTCC"
R1_EXT = "GCGATCGAGGACGGC"
R2_EXT = "CACCGTCTCCGCCTC"
TN5_BC_LEN = 8

R1_PRIMER = R1_EXT + nx.ME
R2_PRIMER = R2_EXT + nx.ME
I1_PRIMER = revcomp(R2_PRIMER)
I2_PRIMER = revcomp(R1_PRIMER)

N5_HEAD = "CACCGCACGAGAGGT"
N5_TAIL = "GTAATCAG"
N7_HEAD = "CAGCACGGCGAGACT"
N7_TAIL = "GACTTGTC"
LIG_BC_LEN = 10
PCR_BC_LEN = 10

SCI3_I1 = nx.ME_RC[-2:] + nx.S7_RC + revcomp(N7_TAIL)
SCI3_I2 = nx.ME_RC[-8:] + nx.S5_RC + revcomp(N5_TAIL)

SOURCE_T5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
SOURCE_2018_P5_1 = "AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGTCGTCGGCAGCGTC"
SOURCE_SCI3_I1 = "CTCCGAGCCCACGAGACGACAAGTC"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def variable(name: str, n: int, letter: str, tag: str, inferred: bool = False) -> Segment:
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred)


def t5(value: str | None = None) -> list[Segment]:
    bc = variable("t5 barcode", TN5_BC_LEN, "A", "cbc") if value is None else seg("t5 barcode", value, "cbc")
    if len(bc) != TN5_BC_LEN: raise ValueError("t5 barcode must be 8 nt")
    return [seg("s5", nx.S5, "s5"), seg("connector A", A_LINK, "r2"), bc,
            seg("Read 1 extension", R1_EXT, "r1"), seg("ME", nx.ME, "me")]


def t7(value: str | None = None) -> list[Segment]:
    bc = variable("t7 barcode", TN5_BC_LEN, "B", "cbc") if value is None else seg("t7 barcode", value, "cbc")
    if len(bc) != TN5_BC_LEN: raise ValueError("t7 barcode must be 8 nt")
    return [seg("s7", nx.S7, "s7"), seg("connector B", B_LINK, "r3"), bc,
            seg("Read 2 extension", R2_EXT, "r1"), seg("ME", nx.ME, "me")]


TWO_LEVEL_PRIMERS = (
    sp.custom("Read 1", "extended mosaic-end A", R1_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
    sp.custom("Index 1 (i7)", "extended index 1", I1_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
    sp.custom("Index 2 (i5)", "extended index 2", I2_PRIMER, "Cusanovich 2018 Table S12"),
    sp.custom("Read 2", "extended mosaic-end B", R2_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
)


def two_level_library(protocol: str) -> Construct:
    if protocol == SCI15:
        pcr_n, inferred = 8, True
    elif protocol == SCI18:
        pcr_n, inferred = 10, False
    else:
        raise ValueError(protocol)
    right = t7()
    con = Construct([
        seg("P5", il.P5, "p5"), variable("PCR i5", pcr_n, "N", "cbc", inferred),
        *t5(), variable("genomic insert", 19, "X", None),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder) for s in reversed(right)],
        variable("PCR i7'", pcr_n, "N", "cbc", inferred), seg("P7'", il.P7_RC, "p7"),
    ], name=f"{protocol} library")
    errors = sp.verify(con, TWO_LEVEL_PRIMERS)
    if errors: raise ValueError("invalid two-level sci-ATAC library: " + "; ".join(errors))
    return con


SCI3_PRIMERS = (
    sp.NEXTERA["R1"],
    sp.custom("Index 1 (i7)", "3LV2 Index 1", SCI3_I1, "Table S7 sequence via secondary source"),
    sp.custom("Index 2 (i5)", "3LV2 Index 2", SCI3_I2, "Table S7 sequence via secondary source"),
    sp.NEXTERA["R2"],
)


def sci3_library() -> Construct:
    """Final sci-ATAC-seq3 library; ligation oligo bases are visibly inferred."""
    inf = True
    con = Construct([
        seg("P5", il.P5, "p5"), variable("i5", PCR_BC_LEN, "I", "cbc"),
        seg("N5 head", N5_HEAD, "r2", inferred=inf),
        variable("N5 barcode", LIG_BC_LEN, "A", "cbc", inf),
        seg("N5 tail", N5_TAIL, "r2", inferred=inf),
        seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"),
        variable("genomic insert", 19, "X", None),
        seg("ME'", nx.ME_RC, "me"), seg("s7'", nx.S7_RC, "s7"),
        seg("N7 tail'", revcomp(N7_TAIL), "r3", inferred=inf),
        variable("N7 barcode'", LIG_BC_LEN, "B", "cbc", inf),
        seg("N7 head'", revcomp(N7_HEAD), "r3", inferred=inf),
        variable("i7'", PCR_BC_LEN, "I", "cbc"), seg("P7'", il.P7_RC, "p7"),
    ], name="sci-ATAC-seq3 library")
    errors = sp.verify(con, SCI3_PRIMERS)
    if errors: raise ValueError("invalid sci-ATAC-seq3 library: " + "; ".join(errors))
    return con


def library(protocol: str) -> Construct:
    return sci3_library() if protocol == SCI3 else two_level_library(protocol)

