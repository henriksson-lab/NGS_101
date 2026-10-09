"""Circular mate-pair junction formation and recovery.

The long molecule is circularised before it is sheared.  Selecting only sheared
fragments that retain the marked circle-closing junction is the defining operation;
ordinary end fragments are deliberately discarded.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, circle_rows


# Illumina's mate-pair technical note prints the common, duplicated junction as
# ME reverse complement followed by ME.  Derive it from the canonical Tn5 ME.
DUPLICATE_JUNCTION = nx.ME_RC + nx.ME


def junction_insert(left: str = "left distant end", right: str = "right distant end",
                    *, duplicate: bool = True) -> list:
    """A recovered shearing fragment spanning the original circle-closing bond."""
    junction = ([seg("junction ME reverse complement", nx.ME_RC, "me"),
                 seg("junction ME", nx.ME, "me")]
                if duplicate else [seg("single junction ME", nx.ME_RC, "me")])
    return [seg(left, "L" * 28, placeholder=True), *junction,
            seg(right, "R" * 28, placeholder=True)]


def recovered_truseq_library() -> tuple[Construct, tuple]:
    """Single-index TruSeq library around a selected mate-pair junction."""
    insert = junction_insert()
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        *insert, seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Nextera mate-pair junction library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("mate-pair library: " + "; ".join(problems))
    return lib, primers


def circularization_rows() -> list[Row]:
    """Long tagged molecule closed at the future mate-pair junction."""
    long = Construct([
        seg("left tag", nx.ME_RC, "me"),
        seg("long genomic fragment", "G" * 48, placeholder=True),
        seg("right tag", nx.ME, "me"),
    ], name="long tagmented genomic molecule")
    return circle_rows(long, "blunt intramolecular ligation")


def tagmentation_scene() -> Scene:
    """Long fragment with both transferred junctions anchored to named segments."""
    parts = [seg("left ME tag", nx.ME_RC, "me"),
             seg("long genomic fragment", "G" * 48, placeholder=True),
             seg("right ME tag", nx.ME, "me")]
    sc = Scene.duplex(parts, label="long tagmented molecule")
    sc.junction("top", "left ME tag", "long genomic fragment", "Tn5 transfer")
    sc.junction("top", "long genomic fragment", "right ME tag", "Tn5 transfer")
    sc.labels("top")
    return sc


def recovery_scene() -> Scene:
    """The short, junction-spanning product retained by streptavidin capture."""
    con = Construct(junction_insert(), name="junction-spanning shear product")
    sc = Scene.duplex(list(con), label="captured shear product")
    sc.junction("top", "junction ME reverse complement", "junction ME",
                "circle-closing junction")
    sc.mark("top", "junction ME reverse complement", "biotin / streptavidin capture")
    sc.labels("top")
    return sc
