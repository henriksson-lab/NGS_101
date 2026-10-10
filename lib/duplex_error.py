"""Physical duplex-preservation models for rare-variant sequencing."""
from __future__ import annotations

from chemdraw import Scene, Segment


def _p(name: str, base: str, n: int, tag: str | None = None) -> Segment:
    return Segment(name, base * n, tag, placeholder=True)


def codec_linked_duplex() -> Scene:
    """One CODEC product containing readable copies of both original strands."""
    parts = [_p("Watson copy", "X", 28), _p("CODEC linkage", "L", 14, "me"),
             _p("Crick copy", "x", 28)]
    sc = Scene(); sc.strand("linked duplex", parts, label="single CODEC cluster molecule")
    sc.junction("linked duplex", "Watson copy", "CODEC linkage", "strand linkage")
    sc.junction("linked duplex", "CODEC linkage", "Crick copy", "strand linkage")
    sc.note("linked duplex", "one read pair carries both complements of the original duplex")
    sc.labels("linked duplex")
    return sc


def nanoseq_restrictive_end_repair() -> Scene:
    """NanoSeq end state: intact restriction ends, with nick extension suppressed."""
    top = [_p("restriction fragment", "X", 34), Segment("blocked nick", "H", placeholder=True)]
    sc = Scene(); sc.strand("top", top, label="restriction-generated genomic duplex")
    sc.anneal("bottom", [_p("restriction fragment partner", "x", 34)], to="top",
              pair=("restriction fragment partner", "restriction fragment"),
              label="restriction-generated genomic duplex")
    sc.mark("top", "blocked nick", "ddNTP terminates nick extension")
    return sc
