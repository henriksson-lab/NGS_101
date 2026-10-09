"""Alternating-pool tiled amplicon schemes.

Neighbouring overlapping amplicons must be placed in different multiplex PCR
reactions.  The constructor enforces that design rule instead of relying on a page
caption or a protocol-specific test.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Row


@dataclass(frozen=True)
class Amplicon:
    name: str
    start: int
    end: int
    pool: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise ValueError(f"{self.name}: invalid interval {self.start}..{self.end}")
        if self.pool not in (1, 2):
            raise ValueError(f"{self.name}: pool must be 1 or 2")


def alternating_scheme(intervals: list[tuple[int, int]]) -> tuple[Amplicon, ...]:
    """Build a two-pool tiling path; adjacent products must overlap."""
    if not intervals:
        raise ValueError("a tiled scheme needs at least one amplicon")
    out = tuple(Amplicon(f"amplicon {i + 1}", a, b, 1 + i % 2)
                for i, (a, b) in enumerate(intervals))
    for left, right in zip(out, out[1:]):
        if right.start <= left.start:
            raise ValueError("amplicons must be ordered by increasing start")
        if right.start >= left.end:
            raise ValueError(f"{left.name} and {right.name} do not overlap")
        if right.pool == left.pool:
            raise ValueError("adjacent overlapping amplicons cannot share a pool")
    return out


def scheme_rows(scheme: tuple[Amplicon, ...], scale: int = 3) -> list[Row]:
    """Render the two multiplex reactions on a common genomic coordinate."""
    if scale < 1:
        raise ValueError("scale must be positive")
    width = max(a.end for a in scheme) // scale + 2
    rows = [Row(chunks=[("viral genome  " + "-" * width, None, False)])]
    for pool in (1, 2):
        chars = [" "] * width
        for a in scheme:
            if a.pool == pool:
                start = a.start // scale
                end = max(start + 1, a.end // scale)
                chars[start] = ">"
                for i in range(start + 1, min(end, width)):
                    chars[i] = "="
                if end < width:
                    chars[end] = ">"
        rows.append(Row(chunks=[(f"PCR pool {pool}   {''.join(chars)}", "cbc", False)]))
    return rows
