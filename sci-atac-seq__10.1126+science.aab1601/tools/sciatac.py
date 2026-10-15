"""Molecular constructors for the individually published sci-ATAC protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

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
SOURCE_SCI3_I1 = "CTCCGAGCCCACGAGACGACAAGTC"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def variable(name: str, n: int, letter: str, tag: str | None,
             inferred: bool = False, **kw) -> Segment:
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred, **kw)


def t5(value: str | None = None) -> list[Segment]:
    semantic = feature("cell_tn5_a", "cell_barcode", "combinatorial",
                       group="cell_id", part="tagmentation A")
    bc = (variable("t5 barcode", TN5_BC_LEN, "A", "cbc", feature=semantic)
          if value is None else seg("t5 barcode", value, "cbc", feature=semantic))
    if len(bc) != TN5_BC_LEN: raise ValueError("t5 barcode must be 8 nt")
    return [seg("s5", nx.S5, "s5"), seg("connector A", A_LINK, "r2"), bc,
            seg("Read 1 extension", R1_EXT, "r1"), seg("ME", nx.ME, "me")]


def t7(value: str | None = None) -> list[Segment]:
    semantic = feature("cell_tn5_b", "cell_barcode", "combinatorial",
                       group="cell_id", part="tagmentation B")
    bc = (variable("t7 barcode", TN5_BC_LEN, "B", "cbc", feature=semantic)
          if value is None else seg("t7 barcode", value, "cbc", feature=semantic))
    if len(bc) != TN5_BC_LEN: raise ValueError("t7 barcode must be 8 nt")
    return [seg("s7", nx.S7, "s7"), seg("connector B", B_LINK, "r3"), bc,
            seg("Read 2 extension", R2_EXT, "r1"), seg("ME", nx.ME, "me")]


TWO_LEVEL_PRIMERS = (
    sp.custom("Read 1", "extended mosaic-end A", R1_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
    sp.custom("Index 1 (i7)", "extended index 1", I1_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
    sp.custom("Index 2 (i5)", "extended index 2", I2_PRIMER, "Cusanovich 2018 Table S12"),
    sp.custom("Read 2", "extended mosaic-end B", R2_PRIMER, "Amini 2014 / Cusanovich 2018 oligo tables"),
)


def two_level_library(protocol: str, *, pcr_n: int,
                      inferred: bool) -> Construct:
    """Assemble either published two-round variant from caller-owned index choices."""
    if pcr_n < 1:
        raise ValueError("PCR index length must be positive")
    right = t7()
    con = Construct([
        seg("P5", il.P5, "p5"), variable("PCR i5", pcr_n, "N", "cbc", inferred,
                                           feature=feature("cell_pcr_i5", "cell_barcode", "combinatorial", group="cell_id", part="PCR i5")),
        *t5(), variable("genomic insert", 19, "X", None),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder, feature=s.feature) for s in reversed(right)],
        variable("PCR i7'", pcr_n, "N", "cbc", inferred,
                 feature=feature("cell_pcr_i7", "cell_barcode", "combinatorial", group="cell_id", part="PCR i7")), seg("P7'", il.P7_RC, "p7"),
    ], name=f"{protocol} library")
    errors = sp.verify(con, TWO_LEVEL_PRIMERS)
    if errors: raise ValueError("invalid two-level sci-ATAC library: " + "; ".join(errors))
    return con


SCI3_PRIMERS = (
    sp.NEXTERA["R1"],
    sp.custom("Index 1 (i7)", "3LV2 Index 1", SCI3_I1, "Domcke 2020 Table S7"),
    sp.custom("Index 2 (i5)", "3LV2 Index 2", SCI3_I2, "Domcke 2020 Table S7"),
    sp.NEXTERA["R2"],
)


def sci3_library() -> Construct:
    """Final sci-ATAC-seq3 library from Domcke 2020 Supplementary Table S7."""
    inf = False
    con = Construct([
        seg("P5", il.P5, "p5"), variable("i5", PCR_BC_LEN, "I", "cbc",
                                          feature=feature("cell_pcr_i5", "cell_barcode", "combinatorial", group="cell_id", part="PCR i5")),
        seg("N5 head", N5_HEAD, "r2", inferred=inf),
        variable("N5 barcode", LIG_BC_LEN, "A", "cbc", inf,
                 feature=feature("cell_ligation_n5", "cell_barcode", "combinatorial", group="cell_id", part="ligation N5")),
        seg("N5 tail", N5_TAIL, "r2", inferred=inf),
        seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"),
        variable("genomic insert", 19, "X", None),
        seg("ME'", nx.ME_RC, "me"), seg("s7'", nx.S7_RC, "s7"),
        seg("N7 tail'", revcomp(N7_TAIL), "r3", inferred=inf),
        variable("N7 barcode'", LIG_BC_LEN, "B", "cbc", inf,
                 feature=feature("cell_ligation_n7", "cell_barcode", "combinatorial", group="cell_id", part="ligation N7")),
        seg("N7 head'", revcomp(N7_HEAD), "r3", inferred=inf),
        variable("i7'", PCR_BC_LEN, "I", "cbc",
                 feature=feature("cell_pcr_i7", "cell_barcode", "combinatorial", group="cell_id", part="PCR i7")), seg("P7'", il.P7_RC, "p7"),
    ], name="sci-ATAC-seq3 library")
    errors = sp.verify(con, SCI3_PRIMERS)
    if errors: raise ValueError("invalid sci-ATAC-seq3 library: " + "; ".join(errors))
    return con


def library(protocol: str) -> Construct:
    if protocol == SCI3:
        return sci3_library()
    if protocol == SCI15:
        return two_level_library(protocol, pcr_n=8, inferred=True)
    raise ValueError(f"{protocol}: supply the paper's two-round index design")


def indexed_transposome(which: str) -> Scene:
    top = t5() if which == "T5" else t7()
    sc = Scene(); sc.strand(which, top, label=which); sc.anneal("pMENTS", [seg("ME'", nx.ME_RC, "me")], to=which, pair=("ME'", "ME"), label="pMENTS", mod5="p")
    sc.mark(which, "t5 barcode" if which == "T5" else "t7 barcode", "round-1 barcode")
    return sc


def tagmented_pair(left: str, right: str) -> Scene:
    """Draw one Tn5 end combination from the same two loaded-end constructors."""
    if left not in {"T5", "T7"} or right not in {"T5", "T7"}:
        raise ValueError("tagmented ends must be T5 or T7")
    g = nx.TAGMENTATION_GAP
    top = [*(t5() if left == "T5" else t7()),
           seg("gap-facing A", "X" * g, placeholder=True),
           seg("core", "XXXXXXXX...XXXXXXXX", placeholder=True)]
    bot = [*(t5() if right == "T5" else t7()),
           seg("gap-facing B", "x" * g, placeholder=True),
           *complement_segments([top[-1]])]
    left_name = f"{left} end" if left != right else f"left {left} end"
    right_name = f"{right} end" if left != right else f"right {right} end"
    sc = Scene(); sc.strand(left_name, top)
    sc.anneal(right_name, bot, to=left_name, pair=("core'", "core"))
    left_ments = f"{left} pMENTS" if left != right else f"left {left} pMENTS"
    right_ments = f"{right} pMENTS" if left != right else f"right {right} pMENTS"
    sc.anneal(left_ments, [seg("ME", nx.ME_RC, "me")],
              to=left_name, pair=("ME", "ME"), label="")
    sc.anneal(right_ments, [seg("ME", nx.ME_RC, "me")],
              to=right_name, pair=("ME", "ME"), label="", above=True)
    sc.mark(left_name, "gap-facing A", "9-nt gap")
    sc.mark(right_name, "gap-facing B", "9-nt gap")
    return sc


def two_tagmented() -> Scene:
    return tagmented_pair("T5", "T7")

def pool_indexed_nuclei(graph: ProtocolGraph, material, *, prefix: str, state, label: str):
    peers = graph.input(prefix + "_peers", label, MaterialKind.MOLECULAR, state=state)
    pooled = graph.product(prefix + "_pool", "pooled " + label,
                           MaterialKind.MOLECULAR, state=state)
    graph.react(prefix + "_pool_step", "Pool indexed nuclei", StepKind.POOL,
                ((material, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    return pooled, [prefix + "_peers", prefix + "_pool_step"]
