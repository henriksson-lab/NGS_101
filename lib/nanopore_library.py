"""Shared topology for motor-adapted Oxford Nanopore libraries.

ONT does not disclose the current adapter bases.  These helpers therefore make the
known molecular roles and junction types structural while preventing callers from
silently presenting proprietary role tokens as nucleotide sequence.
"""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment


def role(name: str, text: str, tag: str | None = None, **kw) -> Segment:
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError("an undisclosed Nanopore component must be a bracketed role token")
    return Segment(name, text, tag, placeholder=True, **kw)


def motor_library(payload: list[Segment], *, adapter: str, name: str) -> Construct:
    """Put the motor-bearing sequencing adapter before a validated payload."""
    if not payload:
        raise ValueError("a Nanopore motor adapter needs a DNA payload")
    return Construct([
        role("motor-loaded sequencing adapter", f"[motor + {adapter}]", "r2"),
        *payload,
    ], name=name)


def single_strand_rows(library: Construct, junctions=()) -> list[Row]:
    """Draw the strand in pore-entry order and mark declared covalent boundaries."""
    sc = Scene()
    sc.strand("pore-entry strand", list(library), label="pore entry")
    for left, right, label in junctions:
        names = [s.name for s in library]
        if left not in names or right not in names or names.index(right) != names.index(left) + 1:
            raise ValueError(f"{left!r} and {right!r} are not an adjacent junction")
        sc.junction("pore-entry strand", left, right, label)
    sc.labels("pore-entry strand")
    return sc.rows()


def no_primer_ending(adapter: str) -> str:
    return (f"No synthesis primer binds. The motor carried by the {adapter} docks at the "
            "nanopore and meters the attached DNA strand through the pore, adapter-bearing "
            "end first.")
