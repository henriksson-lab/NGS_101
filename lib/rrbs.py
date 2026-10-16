"""Shared RRBS model with explicitly inferred canonical adapter arms.

The scRRBS article names standard premethylated indexed Illumina adapters but
does not print their exact oligos or selected index member.  Canonical TruSeq
arms are therefore a conditional illustration in each inheriting protocol.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail


INDEX_NT = 6
MSPI_SITE = "CCGG"
MSPI_OVERHANG = "CG"
CELL_INDEX = feature("cell_i7", "cell_barcode", "whitelist",
                     whitelist="published indexed primer set")
UNIVERSAL = il.TRUSEQ_P5_FULL
INDEXED_FIXED_5 = il.INDEX1_PRIMER


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def universal_segments(*, inferred: bool = True) -> list[Segment]:
    return [
        _seg("P5", il.P5, "p5", inferred=inferred),
        _seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1", inferred=inferred),
    ]


def indexed_segments(*, inferred: bool = True) -> list[Segment]:
    return [
        _seg("Index-primer / Read 2 arm", INDEXED_FIXED_5, "r2", inferred=inferred),
        _seg("i7 cell index", "N" * INDEX_NT, "cbc", placeholder=True,
             inferred=inferred, feature=CELL_INDEX),
        _seg("P7'", il.P7_RC, "p7", inferred=inferred),
    ]


def adapter_scene() -> Scene:
    stem_n = len(il.STEM)
    universal = [
        _seg("universal fork", UNIVERSAL[:-(stem_n + 1)], "p5", inferred=True),
        _seg("universal stem", UNIVERSAL[-(stem_n + 1):-1], "r1", inferred=True),
        _seg("3' T", UNIVERSAL[-1], "r1", inferred=True),
    ]
    indexed = [
        _seg("indexed stem", INDEXED_FIXED_5[:stem_n], "r2", inferred=True),
        _seg("indexed fork", INDEXED_FIXED_5[stem_n:], "r2", inferred=True),
        _seg("i7", "N" * INDEX_NT, "cbc", placeholder=True, inferred=True,
             feature=CELL_INDEX),
        _seg("P7'", il.P7_RC, "p7", inferred=True),
    ]
    sc = Scene()
    sc.strand("universal", universal, label="universal")
    sc.anneal("indexed", indexed, to="universal",
              pair=("indexed stem", "universal stem"), label="indexed", mod5="p",
              unpaired=("indexed fork", "i7", "P7'"))
    sc.mark("universal", "universal stem", "12-bp TruSeq stem")
    return sc


def digested_fragment_scene() -> Scene:
    core = [_seg("paired fragment", "GXXXXXXXX...XXXXXXXXCC", placeholder=True)]
    sc = Scene()
    sc.strand("top", [_seg("left 5' CG", MSPI_OVERHANG), *core], label="fragment")
    sc.anneal("bottom", [_seg("right 5' CG", MSPI_OVERHANG),
                         *complement_segments(core)],
              to="top", pair=("paired fragment'", "paired fragment"), label="fragment",
              mod5="p")
    sc.mark("top", "left 5' CG", "MspI 5' overhang")
    sc.mark("bottom", "right 5' CG", "MspI 5' overhang")
    return sc


def a_tailed_fragment_scene() -> Scene:
    core = Construct([_seg("filled fragment", "CGGXXXXXXXX...XXXXXXXXCCG",
                           placeholder=True)])
    return dA_tailed_scene(repair_and_dA_tail(core), label="fragment")


def final_library() -> Construct:
    return Construct([
        *universal_segments(),
        _seg("bisulfite-converted insert", "YGGXXXXXXXX...XXXXXXXXYTG",
             placeholder=True),
        _seg("dA junction", "A"),
        *indexed_segments(),
    ], name="scRRBS final library")


SEQ_PRIMERS = [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"]]


def _validate_model() -> None:
    if UNIVERSAL != il.P5 + il.TRUSEQ_READ1[4:]:
        raise ValueError("universal adaptor no longer assembles from canonical P5/Read 1")
    if revcomp(il.TRUSEQ_READ2) != "A" + INDEXED_FIXED_5:
        raise ValueError("A-tail plus indexed adaptor no longer forms the Read 2 site")
    adapter_scene().rows()
    problems = sp.verify(final_library(), SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scRRBS sequencing layout: " + "; ".join(problems))


_validate_model()
