"""Molecular constructors for the published inDrop v1 and v2 protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import seqprimers as sp

V1 = "inDrop v1"
V2 = "inDrop v2"

LEADER = "CGATGACG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
MID = "ATACCACCATGG"
PE1 = il.TRUSEQ_READ1[3:]
W1_STAR = "AAGGCGTCACAAGCAATCACTC"
W1 = revcomp(W1_STAR)
BC2_LEN = 8
UMI_LEN = 6
DT_LEN = 19
PE2 = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCT"

RLO = revcomp(PE2[2:])
V1_R2 = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
V1_R1 = il.TRUSEQ_READ1[4:]
V2_R1 = PE2[2:]
V2_INDEX = revcomp(PE1)
V2_R2 = PE1

SOURCE_ACRYDITE_BASES = "CGATGACGTAATACGACTCACTATAGGGATACCACCATGGCTCTTTCCCTACACGACGCTCTTC"
SOURCE_PE2_N6 = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNNN"
BC1_FEATURE = feature("cell_barcode_1", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 1")
BC2_FEATURE = feature("cell_barcode_2", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 2")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"barcode 1": BC1_FEATURE, "barcode 1'": BC1_FEATURE,
            "barcode 2": BC2_FEATURE, "barcode 2'": BC2_FEATURE,
            "UMI": UMI_FEATURE, "i7 read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def var(name: str, n: int, letter: str, tag: str | None = None,
        inferred: bool = False) -> Segment:
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred)


def bead_primer(bc1_len: int = 8, *, inferred: bool = False) -> Construct:
    if bc1_len not in (8, 9, 10, 11):
        raise ValueError("inDrop barcode 1 must be 8-11 nt")
    return Construct([
        seg("leader", LEADER, inferred=inferred), seg("T7 promoter", T7_PROMOTER, "tso", inferred=inferred),
        seg("middle", MID, inferred=inferred), seg("PE1", PE1, "r2", inferred=inferred),
        var("barcode 1'", bc1_len, "A", "cbc", inferred), seg("W1", W1, "r3", inferred=inferred),
        var("barcode 2'", BC2_LEN, "B", "cbc", inferred), var("UMI", UMI_LEN, "N", "umi", inferred),
        seg("T19", "T" * DT_LEN, inferred=inferred), var("V", 1, "V", None, inferred),
    ], name="finished inDrop bead primer")


V1_PRIMERS = (
    sp.custom("Read 1", "v1 Read 1", V1_R1, "upstream secondary source"),
    sp.custom("Read 2", "v1 Read 2", V1_R2, "upstream secondary source"),
)
V2_PRIMERS = (
    sp.custom("Read 1", "Custom Read 1", V2_R1, "Nature Protocols supplementary oligo table"),
    sp.custom("Index 1 (i7)", "Custom Index Read", V2_INDEX, "Nature Protocols supplementary oligo table"),
    sp.custom("Read 2", "Custom Read 2", V2_R2, "Nature Protocols supplementary oligo table"),
)


def v1_library(bc1_len: int = 8) -> Construct:
    """Entire v1 sequence model is secondary-source and therefore inseparably inferred."""
    inf = True
    con = Construct([
        seg("P5", il.P5, "p5", inferred=inf),
        seg("Read 1 site", V1_R1, "r1", inferred=inf),
        var("barcode 1'", bc1_len, "A", "cbc", inf), seg("W1", W1, "r3", inferred=inf),
        var("barcode 2'", BC2_LEN, "B", "cbc", inf), var("UMI", UMI_LEN, "N", "umi", inf),
        seg("T19", "T" * DT_LEN, inferred=inf), var("V", 1, "V", None, inf),
        var("antisense mRNA", 19, "X", "r1", inf), seg("RLO", RLO, "r2", inferred=inf),
        seg("RT2 extension", "GAGACCG", "r2", inferred=inf), seg("P7'", il.P7_RC, "p7", inferred=inf),
    ], name="inDrop v1 library")
    _validate(con, V1_PRIMERS)
    return con


def v2_library(bc1_len: int = 8) -> Construct:
    if bc1_len not in (8, 9, 10, 11): raise ValueError("barcode 1 must be 8-11 nt")
    con = Construct([
        seg("P5", il.P5, "p5"), seg("GGTC", "GGTC", "r1"), seg("PE2", PE2, "r1"),
        var("sense mRNA", 19, "X", "r1"), seg("A19", "A" * DT_LEN), var("B", 1, "B"),
        var("UMI", UMI_LEN, "N", "umi"), var("barcode 2", BC2_LEN, "B", "cbc"),
        seg("W1*", W1_STAR, "r3"), var("barcode 1", bc1_len, "A", "cbc"),
        seg("PE1*", revcomp(PE1), "r2"), var("i7 read", 6, "I", "cbc"),
        seg("P7'", il.P7_RC, "p7"),
    ], name="inDrop v2 library")
    _validate(con, V2_PRIMERS)
    return con


def _validate(con: Construct, primers) -> None:
    for p in primers:
        if sp.locate(con, p) is None:
            raise ValueError(f"{con.name}: {p.name} site missing")


def library(protocol: str) -> Construct:
    if protocol == V1: return v1_library()
    if protocol == V2: return v2_library()
    raise ValueError(protocol)

