"""Small supplied-state drawing helpers used by protocol-owned models."""
from __future__ import annotations

from chemdraw import Scene, Segment, feature


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


# Pass as ``nt`` when the authoritative source does not state an index
# length: the window is drawn as a bracketed role token, which carries no molecular
# length, instead of an assumed run of index letters.
UNASSIGNED_LENGTH = None


def index_window(name: str, letter: str, nt: int | None, identifier: str, *,
                 inferred: bool = False, sample_feature=None) -> Segment:
    """A sample-index window of ``nt`` placeholder bases, or of unassigned length."""
    if nt is UNASSIGNED_LENGTH:
        return seg(name, f"[{name.split()[0]} index; length not stated]", "cbc",
                   placeholder=True,
                   inferred=inferred,
                   feature=sample_feature or feature(identifier, "sample_index", "unknown",
                                                     note="index length not stated by the source"))
    return seg(name, letter * nt, "cbc", placeholder=True, inferred=inferred,
               feature=sample_feature or feature(identifier, "sample_index", "unknown"))


def duplex_fragment(label: str = "genomic fragment", length: int = 34) -> Scene:
    return Scene.duplex([seg(label, "X" * length, placeholder=True)], label=label)


def joined_scene(parts: list[Segment], junctions: tuple[tuple[str, str, str], ...],
                 *, label: str, duplex: bool = True) -> Scene:
    """Draw covalently joined parts; ``Scene.junction`` validates every boundary."""
    if duplex:
        sc = Scene.duplex(parts, label=label)
        strand = "top"
    else:
        sc = Scene()
        strand = "molecule"
        sc.strand(strand, parts, label=label)
    for left, right, kind in junctions:
        sc.junction(strand, left, right, kind)
    sc.labels(strand)
    return sc
