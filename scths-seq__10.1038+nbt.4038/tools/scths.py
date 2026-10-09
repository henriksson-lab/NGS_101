"""Molecular construct model for scTHS-seq (Lake et al., 2018)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


T7_LEADER = "AATTAATACGACTCACTATA"
CONNECTOR = "GGGAGATCCACGCGC"
R5_EXAMPLE = "TCTAAT"
I5_OLIGO_INDEX = "CTCTCTAT"
I7_OLIGO_INDEX = "TCGCCTTA"
I7_INDEX_READ = "TAAGGCGA"


def r5_transposon() -> list[Segment]:
    return [seg("T7 leader/promoter", T7_LEADER, "t7"),
            seg("IVT/i5 connector", CONNECTOR, "t7"),
            seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def s7_transposon() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]


def first_tagmented_scene(insert_nt: int = 28) -> Scene:
    """One end of a genomic fragment tagged by the T7/r5/s5 homodimer."""
    top = [*r5_transposon(), seg("accessible genomic DNA", "X" * insert_nt,
                                placeholder=True)]
    sc = Scene()
    sc.strand("top", top, label="round-1 tagged DNA")
    sc.anneal("ME bottom", [seg("mosaic end'", nx.ME_RC, "me")], to="top",
              pair=("mosaic end'", "mosaic end"), label="", mod5="p")
    sc.mark("top", "mosaic end", "9-nt genomic gap on the opposite strand")
    return sc


def ivt_template(insert_nt: int = 28) -> Construct:
    return Construct([seg("IVT/i5 connector", CONNECTOR, "t7"),
                      seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True),
                      seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me"),
                      seg("genomic DNA", "X" * insert_nt, placeholder=True)],
                     name="scTHS-seq re-made dsDNA")


INDEX2_P5 = sp.custom(
    "Index 2 (i5)", "scTHS nXTv2 i5 index-read primer", il.P5,
    "Lake et al. 2018 Supplementary Table S12",
    "Reads i5, the 15-nt connector, the r5 barcode and three s5 bases.")
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], INDEX2_P5)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 index", I5_OLIGO_INDEX, "cbc",
                                        feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        seg("IVT/i5 connector", CONNECTOR, "t7"),
        seg("round-1 r5 barcode", "B" * 6, "cbc", placeholder=True,
            feature=feature("cell_r5", "cell_barcode", "combinatorial", group="cell_id", part="round 1 r5")),
        seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("second mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc",
            feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scTHS-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final scTHS-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")
