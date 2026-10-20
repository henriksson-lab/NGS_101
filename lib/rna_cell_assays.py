"""Reusable state, scene, and library primitives for RNA and cell assays."""
from __future__ import annotations


from chemdraw import Scene, Segment


def single(label: str, parts: tuple[Segment, ...], *, note: str = "",
           junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    """One selectable strand with checked, named covalent junctions."""
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, caption in junctions:
        sc.junction("molecule", left, right, caption)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc
