"""Molecular construct model for scifi-RNA-seq (Datlinger et al.)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp

UMI_FEATURE = feature("umi", "umi", "random")
ROUND1_FEATURE = feature("cell_barcode_round1", "cell_barcode", "combinatorial",
                         group="cell_barcode", part="RT plate")
ROUND2_FEATURE = feature("cell_barcode_round2", "cell_barcode", "combinatorial",
                         group="cell_barcode", part="bead")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"UMI": UMI_FEATURE, "round-1 barcode": ROUND1_FEATURE,
            "round-2 bead barcode": ROUND2_FEATURE, "i7 index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


RT_BARCODE_EXAMPLE = "AGTGATTAGCA"
P7_INDEX_READ = "TAAGGCGA"
BRIDGE = "CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA"
PARTIAL_P5 = il.P5[:22]
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"


def rt_primer() -> list[Segment]:
    return [seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
            seg("UMI", "N" * 8, "umi", placeholder=True), seg("fixed N", "A"),
            seg("round-1 barcode", "B" * 11, "cbc", placeholder=True),
            seg("fixed V", "A"), seg("poly(T)", "T" * 30),
            seg("V", "V", placeholder=True), seg("N", "N", placeholder=True)]


def bead_oligo() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("round-2 bead barcode", "C" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def bridge_oligo() -> list[Segment]:
    return [seg("RT-primer bridge", BRIDGE[:21]), seg("bead bridge", BRIDGE[21:])]


def ligation_scene() -> Scene:
    """Bridge-enforced nick ligation between bead and phosphorylated RT oligos."""
    primer = rt_primer()
    joined = [*bead_oligo(),
              seg("Read-1 bridge site", il.TRUSEQ_READ1[:21], "r1"),
              seg("Read-1 remainder", il.TRUSEQ_READ1[21:], "r1"), *primer[1:]]
    sc = Scene()
    sc.strand("product", joined, label="bead + RT primer")
    sc.anneal("bridge", bridge_oligo(), to="product",
              pair=("RT-primer bridge", "Read-1 bridge site"), label="3'-ddC bridge",
              above=True)
    # The named pair fixes the 21-nt register; the adjacent 14 nt necessarily cover s5
    # because BRIDGE is the reverse complement of s5 + the first 21 nt of Read 1.
    if revcomp(BRIDGE) != nx.S5 + il.TRUSEQ_READ1[:21]:
        raise ValueError("bridge no longer joins s5 directly to TruSeq Read 1")
    sc.junction("product", "s5", "Read-1 bridge site")
    return sc


def ligated_first_strand(insert_nt: int = 28) -> Construct:
    return Construct([
        *bead_oligo(), *rt_primer(),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("untemplated C", "CCC"),
    ], name="scifi-RNA ligated first strand")


def tso_scene(insert_nt: int = 28) -> Scene:
    first = ligated_first_strand(insert_nt)
    tso = [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
           seg("rGrGrG", "rGrGrG", placeholder=True)]
    sc = Scene()
    sc.strand("first", list(first), label="first-strand cDNA")
    sc.anneal("TSO", tso, to="first", pair=("rGrGrG", "untemplated C"), label="TSO")
    sc.arrow("TSO", "template switching; then TSO-primer enrichment")
    return sc


INDEX2_FROM_P5 = sp.custom(
    "Index 2 (i5)", "flow-cell P5 oligo (forward-strand workflow)", il.P5,
    "Illumina forward-strand indexed sequencing",
    "The bead barcode occupies the i5 position immediately after P5.")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], INDEX2_FROM_P5, sp.NEXTERA["R2"])


def final_library(insert_nt: int = 28) -> Construct:
    """Selected bead-end/i7-Tn5 product after partial-P5/P7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("round-2 bead barcode", "C" * 16, "cbc", placeholder=True),
        seg("s5", nx.S5, "s5"), seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("fixed N", "A"),
        seg("round-1 barcode", "B" * 11, "cbc", placeholder=True),
        seg("fixed V", "A"), seg("poly(T)", "T" * 30),
        seg("V", "V", placeholder=True), seg("N", "N", placeholder=True),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", P7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scifi-RNA-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final scifi-RNA-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":21})
    by_id={x.feature.id:(x.cycles,x.feature.label) for x in spans}
    return [by_id["umi"],("9","fixed N"),by_id["cell_barcode_round1"],
            ("21","fixed V")]
