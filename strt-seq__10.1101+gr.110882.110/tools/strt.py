"""Constructors for the three individually published STRT protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, revcomp
import illumina as il
import nextera as nx
import rt
import seqprimers as sp

STRT = "STRT-seq"
C1 = "STRT-C1"
TWO_I = "STRT-seq-2i"
PROTOCOLS = (STRT, C1, TWO_I)

H = rt.SMART_HANDLE
STRT_DT_LINK = "CGAC"
STRT_TSO_LINK = "GCAGTGCT"
STRT_BC_LEN = 6
STRT_DT_LEN = 30

C1_HANDLE = il.P5[:20] + "T"
C1_DT_LINK = "CG"
C1_DT_LEN = 31
C1_UMI_LEN = 5
C1_PCR = "G" + C1_HANDLE

P1B = il.TRUSEQ_READ1[-22:]
TWO_I_UMI_LEN = 6
WELL_INDEX_LEN = 5
SUBARRAY_INDEX_LEN = 8

TN5_TAIL = nx.S5[-5:] + nx.ME
TN5_U = revcomp(TN5_TAIL)

STRT_DT = H + STRT_DT_LINK + "T" * STRT_DT_LEN + "VN"
STRT_TSO = H + STRT_TSO_LINK + "N" * STRT_BC_LEN + "GGG"
P2_TOP = il.P7[:22] + il.TRUSEQ_READ2[-12:]
STRT_LIB_PCR1 = il.P5[:25] + H
STRT_LIB_PCR2 = il.P7[:22]

C1_DT = C1_HANDLE + C1_DT_LINK + "T" * C1_DT_LEN
C1_TSO = C1_HANDLE + "N" * C1_UMI_LEN + "GGG"
C1_TN5 = il.P7[:21] + "N" * SUBARRAY_INDEX_LEN + TN5_TAIL

TWO_I_TSO = P1B + "N" * TWO_I_UMI_LEN + "GGG"
TWO_I_WELL_PRIMER = il.P5 + "N" * WELL_INDEX_LEN + P1B[:-1]
TWO_I_READ1 = il.P5[1:] + "N" * 6 + P1B
TWO_I_INDEX2 = il.P5

# Verbatim source rows used only at the source/model boundary checks.
SOURCE_STRT_V3 = "AAGCAGTGGTATCAACGCAGAGTCGACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN"
SOURCE_C1_TN5_1 = "CAAGCAGAAGACGGCATACGACGTCTAATGCGTCAGATGTGTATAAGAGACAG"
SOURCE_TWO_I_WELL = "AATGATACGGCGACCACCGAGATCTACACXXXXXCTACACGACGCTCTTCCGATC"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def variable(name: str, n: int, letter: str, tag: str, *, inferred: bool = False) -> Segment:
    if n < 1:
        raise ValueError(f"{name} must have positive length")
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred)


def tn5_indexed_top(value: str | None = None, *, inferred: bool = False) -> list[Segment]:
    bc = variable("subarray barcode", SUBARRAY_INDEX_LEN, "S", "cbc", inferred=inferred)
    if value is not None:
        if len(value) != SUBARRAY_INDEX_LEN:
            raise ValueError("subarray barcode must be 8 nt")
        bc = seg("subarray barcode", value, "cbc", inferred=inferred)
    return [seg("P7 fragment", il.P7[:21], "p7", inferred=inferred),
            bc,
            seg("GCGTC", nx.S5[-5:], "s5", inferred=inferred),
            seg("ME", nx.ME, "me", inferred=inferred)]


def strt_library() -> Construct:
    con = Construct([
        seg("P5 fragment", il.P5[:25], "p5"), seg("SMART handle", H, "tso"),
        seg("TSO spacer", STRT_TSO_LINK, "r1"),
        variable("cell barcode", STRT_BC_LEN, "B", "cbc"), seg("GGG", "GGG", "tso"),
        variable("mRNA 5' end", 19, "X", "r1"),
        seg("P2 adapter'", revcomp(P2_TOP), "p7"),
    ], name="STRT-seq library")
    if sp.locate(con, STRT_READ1) is None:
        raise ValueError("STRT custom Read 1 site is absent")
    return con


def c1_library() -> Construct:
    right = tn5_indexed_top()
    con = Construct([
        seg("extra G", "G", "p5"), seg("C1 handle", C1_HANDLE, "p5"),
        variable("UMI", C1_UMI_LEN, "N", "umi"), seg("GGG", "GGG", "tso"),
        variable("mRNA 5' end", 19, "X", "r1"),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder) for s in reversed(right)],
    ], name="STRT-C1 library")
    for primer in C1_PRIMERS:
        if sp.locate(con, primer) is None:
            raise ValueError(f"{primer.name} site is absent")
    return con


def two_i_library() -> Construct:
    """One explicit outcome of the unpublished P2-primer/barcode junction.

    The full subarray barcode followed by P7' is retained, but both are styled inferred
    because the P2 primer's three terminal bases overlap the barcode and the paper does
    not resolve whether proofreading trims or overwrites them.
    """
    right = tn5_indexed_top(inferred=True)
    con = Construct([
        seg("P5", il.P5, "p5"), variable("well index", WELL_INDEX_LEN, "W", "cbc"),
        seg("P1B", P1B, "r1"), variable("UMI", TWO_I_UMI_LEN, "N", "umi"),
        seg("GGG", "GGG", "tso"), variable("mRNA 5' end", 19, "X", "r1"),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder, True) for s in reversed(right)],
        seg("P7'", il.P7_RC, "p7", inferred=True),
    ], name="STRT-seq-2i library")
    for primer in TWO_I_INDEX_PRIMERS:
        if sp.locate(con, primer) is None:
            raise ValueError(f"{primer.name} site is absent")
    return con


STRT_READ1 = sp.custom("Read 1", "derived STRT custom-primer site",
                       il.P5[20:25] + H + STRT_TSO_LINK,
                       "The paper states custom primer but does not print it; site reconstructed from the library")
C1_READ1 = sp.custom("Read 1", "C1 handle site", C1_HANDLE,
                     "Landing site derived from the published construct; primer sequence not printed")
C1_INDEX1 = sp.custom("Index 1 (i7)", "C1-TN5-U", TN5_U,
                      "Supplementary oligo table")
C1_PRIMERS = (C1_READ1, C1_INDEX1)
TWO_I_INDEX1 = sp.custom("Index 1 (i7)", "STRT-TN5-U", TN5_U,
                         "Hochgerner et al. Table 1")
TWO_I_INDEX2_PRIMER = sp.custom("Index 2 (i5)", "DI_idxP1A-Seq", TWO_I_INDEX2,
                                "Hochgerner et al. Table 1")
TWO_I_INDEX_PRIMERS = (TWO_I_INDEX1, TWO_I_INDEX2_PRIMER)


def library(protocol: str) -> Construct:
    if protocol == STRT:
        return strt_library()
    if protocol == C1:
        return c1_library()
    if protocol == TWO_I:
        return two_i_library()
    raise ValueError(f"unknown STRT protocol {protocol!r}")
