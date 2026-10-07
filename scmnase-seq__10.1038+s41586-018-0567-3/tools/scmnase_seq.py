"""Molecular construct model for scMNase-seq.

The paper specifies universal adapters and indexed primers but prints no oligo bases.
Every adapter-derived segment is therefore constructed through :func:`inferred`.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments

INDEX_NT = 6


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """An adapter region absent from the defining paper, always dotted."""
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


def digested_fragment_scene() -> Scene:
    """Representative nucleosome-protected fragment left by MNase."""
    return Scene.duplex([seg("MNase-protected genomic fragment", "X" * 44,
                             placeholder=True)], label="protected DNA")


def a_tailed_fragment_scene() -> Scene:
    """End-repaired fragment with one unpaired 3'-A at each end."""
    core = [inferred("end-repaired fragment", "X" * 40, placeholder=True)]
    sc = Scene()
    sc.strand("top", [*core, inferred("right dA", "A")], label="fragment")
    sc.anneal("bottom", [*complement_segments(core), inferred("left dA", "A")],
              to="top", pair=("end-repaired fragment'", "end-repaired fragment"),
              label="fragment", unpaired=("left dA",))
    sc.mark("top", "right dA", "3'-A")
    sc.mark("bottom", "left dA", "3'-A")
    return sc


def adapter_bottom() -> list[Segment]:
    return [
        inferred("Read 1 arm", il.TRUSEQ_READ1[:-(len(il.STEM_COMPLEMENT) + 1)], "r1"),
        inferred("stem complement", il.STEM_COMPLEMENT, "r1"),
        inferred("3' T", "T", "r1"),
    ]


def adapter_top() -> list[Segment]:
    return [
        inferred("stem", il.STEM, "r2"),
        inferred("unpaired Read 2 arm", il.INDEX1_PRIMER[len(il.STEM):20], "r2"),
    ]


def adapter_scene() -> Scene:
    """Inferred TruSeq fork; Scene enforces the 12-bp stem."""
    sc = Scene()
    sc.strand("bottom", adapter_bottom(), label="3'-T strand")
    sc.anneal("top", adapter_top(), to="bottom", pair=("stem", "stem complement"),
              label="5'-phosphorylated strand", mod5="p",
              unpaired=("unpaired Read 2 arm",))
    sc.mark("bottom", "3' T", "ligation overhang")
    return sc


def p5_primer() -> list[Segment]:
    return [inferred("P5", il.P5, "p5"),
            inferred("Read 1 continuation", il.TRUSEQ_READ1[4:], "r1")]


def indexed_p7_primer() -> list[Segment]:
    return [inferred("P7", il.P7, "p7"),
            inferred("i7", "N" * INDEX_NT, "cbc", placeholder=True),
            inferred("Read 2", il.TRUSEQ_READ2, "r2")]


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 44) -> Construct:
    """PCR-completed single-index library; required primer sites are validated here."""
    lib = Construct([
        inferred("P5", il.P5, "p5"),
        inferred("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("MNase-protected genomic insert", "X" * insert_nt, placeholder=True),
        inferred("dA junction", "A"),
        inferred("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        inferred("i7 reverse complement", "N" * INDEX_NT, "cbc", placeholder=True),
        inferred("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scMNase-seq library")
    for primer in SEQ_PRIMERS:
        if sp.locate(lib, primer) is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
    return lib
