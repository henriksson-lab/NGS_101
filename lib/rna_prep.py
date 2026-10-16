"""Reusable molecular states for RNA-first library construction."""
from __future__ import annotations

from chemdraw import Scene, Segment


def seg(name: str, sequence: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, sequence, tag, **kw)


def rna_template(name: str = "RNA insert", length: int = 34) -> list:
    """A typed RNA starting state; X is a length-preserving unknown sequence."""
    sc = Scene()
    sc.strand("RNA", [seg(name, "X" * length, placeholder=True)], label="RNA 5'→3'")
    sc.labels("RNA")
    return sc.rows()


def rna_dna_hybrid(name: str = "RNA/cDNA hybrid", length: int = 34) -> list:
    """An RNA/cDNA heteroduplex with pairing enforced by :class:`Scene`."""
    top = [seg("RNA template", "X" * length, placeholder=True)]
    return Scene.duplex(top, label=name).rows()


def single_strand(name: str, segments: list[Segment], label: str = "") -> list:
    sc = Scene()
    sc.strand(name, segments, label=label or name)
    sc.labels(name)
    return sc.rows()
