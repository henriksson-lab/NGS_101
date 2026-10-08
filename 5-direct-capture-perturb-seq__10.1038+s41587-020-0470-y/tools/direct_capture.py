"""5-prime direct-capture Perturb-seq (standard sgRNA-CR1 branch)."""
from __future__ import annotations

import crispr
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


TENX_RT_ADAPTER = "AAGCAGTGGTATCAACGCAGAGTAC"
GUIDE_ANNEAL = "CAAGTTGATAACGGACTAGCC"
OJR160 = TENX_RT_ADAPTER + GUIDE_ANNEAL

OJR163 = il.TRUSEQ_P5_FULL
I7 = "AGGAGTCC"                         # example index printed in the protocol
NESTED_ANNEAL = TENX_RT_ADAPTER[-5:] + GUIDE_ANNEAL
OJR165 = il.P7 + I7 + nx.ADAPTOR_S7 + NESTED_ANNEAL

SCAFFOLD = crispr.SCAFFOLD_FE             # the paper's optimized CR1 constant region
SCAFFOLD_SITE = revcomp(GUIDE_ANNEAL)
SITE_OFFSET = SCAFFOLD.index(SCAFFOLD_SITE)

PARTIAL_R1 = il.TRUSEQ_READ1[11:]
SWITCH_SPACER = "TTTCTTATAT"


def ojr160_segments():
    return [seg("10x RT adapter", TENX_RT_ADAPTER, "r1"),
            seg("guide annealing", GUIDE_ANNEAL, "tso")]


def ojr163_segments():
    return [seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1")]


def ojr165_segments():
    return [seg("P7", il.P7, "p7"), seg("i7", I7, "cbc"),
            seg("Read 2 arm", nx.ADAPTOR_S7, "r2"),
            seg("nested guide site", NESTED_ANNEAL, "tso")]


def bead_tso():
    return [
        seg("Partial Read 1", PARTIAL_R1, "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True),
        seg("UMI", "U" * 10, "umi", placeholder=True),
        seg("switch spacer", SWITCH_SPACER, "tso"),
        seg("rGrGrG", "GGG", "tso"),
    ]


def guide_rna():
    return [
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 before capture site", SCAFFOLD[:SITE_OFFSET], "insert"),
        seg("oJR160 site", SCAFFOLD_SITE, "tso"),
        seg("CR1 remainder", SCAFFOLD[SITE_OFFSET + len(SCAFFOLD_SITE):], "insert"),
    ]


def capture_scene():
    sc = Scene()
    sc.strand("guide", guide_rna(), label="sgRNA-CR1")
    sc.anneal(
        "oJR160",
        ojr160_segments(),
        to="guide", pair=("guide annealing", "oJR160 site"),
        label="oJR160", unpaired=("10x RT adapter",),
    )
    sc.mark("guide", "protospacer", "guide identity")
    sc.mark("guide", "oJR160 site", "guide-specific RT site")
    sc.arrow("oJR160", "reverse transcriptase")
    return sc


def first_strand():
    return Construct([
        *ojr160_segments(),
        seg("copied CR1", revcomp(SCAFFOLD[:SITE_OFFSET]), "insert"),
        seg("copied protospacer", "n" * 20, "insert", placeholder=True),
        seg("CCC", "CCC", "tso"),
    ], name="guide first strand")


def switch_scene():
    first = first_strand()
    sc = Scene()
    sc.strand("first", list(first), label="guide first strand")
    sc.anneal("bead TSO", bead_tso(), to="first", pair=("rGrGrG", "CCC"),
              label="barcoded bead TSO", above=True, mod5="bead")
    sc.mark("bead TSO", "cell barcode", "cell barcode + UMI", through="UMI")
    sc.arrow("first", "RT switches template")
    return sc


def amplified_guide_cdna():
    """Guide-bearing 10x cDNA before the guide-specific library PCR."""
    return Construct([
        *bead_tso(),
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 through capture site", SCAFFOLD[:SITE_OFFSET + len(SCAFFOLD_SITE)],
            "insert"),
        seg("RT adapter complement", revcomp(TENX_RT_ADAPTER), "r1"),
    ], name="barcoded guide cDNA")


def final_library():
    """The approximately 250-bp guide library produced by oJR163/oJR165 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True),
        seg("UMI", "U" * 10, "umi", placeholder=True),
        seg("switch spacer", SWITCH_SPACER, "tso"),
        seg("GGG", "GGG", "tso"),
        seg("protospacer", "N" * 20, "insert", placeholder=True),
        seg("CR1 before capture site", SCAFFOLD[:SITE_OFFSET], "insert"),
        seg("oJR160 site", SCAFFOLD_SITE, "tso"),
        seg("nested adapter site", revcomp(TENX_RT_ADAPTER[-5:]), "r1"),
        seg("Index 1 / Read 2 arm", nx.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", revcomp(I7), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="5-prime direct-capture guide library")
    problems = sp.verify(lib, sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid direct-capture library: " + "; ".join(problems))
    return lib


def final_scene():
    lib = final_library()
    sc = Scene.duplex(list(lib), label="guide library")
    sc.mark("top", "cell barcode", "Read 1: cell barcode", through="UMI")
    sc.mark("top", "protospacer", "guide identity reached by Read 2")
    sc.labels("top")
    return sc


def sequencing_primers():
    return (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def read_layout():
    return [
        ("Read 1", "1–16", "cell barcode"),
        ("Read 1", "17–26", "UMI"),
        ("Index 1", "1–8", "sample index"),
        ("Read 2", "98 cycles recommended", "constant region, then protospacer"),
    ]
