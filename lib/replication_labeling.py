"""Shared molecular states for bulk and single-cell replication-timing assays."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment


def _part(name: str, width: int = 20) -> Segment:
    return Segment(name, "X" * width, placeholder=True)


@dataclass(frozen=True)
class NascentPulse:
    analogue: str
    fractions: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.analogue:
            raise ValueError("nascent-DNA pulse requires a named nucleotide analogue")
        if len(self.fractions) < 2:
            raise ValueError("replication timing requires at least two temporal fractions")

    def labelled_scene(self) -> Scene:
        sc = Scene(); sc.strand("nascent", [_part("parental DNA"),
                                             _part(f"nascent DNA containing {self.analogue}")],
                                      label="pulse-labelled replicating DNA")
        sc.note("nascent", f"{self.analogue} is incorporated only during the pulse")
        sc.labels("nascent")
        return sc

    def sorted_scene(self) -> Scene:
        sc = Scene()
        for i, fraction in enumerate(self.fractions):
            sc.strand(f"fraction {i}", [_part(f"{fraction} nascent DNA")], label=fraction)
        sc.note("fraction 0", "fractions are sorted by DNA content before BrdU capture")
        return sc


@dataclass(frozen=True)
class CopyNumberReplication:
    cell_fraction: str
    amplification: str

    def __post_init__(self) -> None:
        if "S" not in self.cell_fraction:
            raise ValueError("single-cell replication profiling requires an S-phase cell")
        if not self.amplification:
            raise ValueError("single-cell copy-number profiling requires WGA")

    def states(self) -> Scene:
        sc = Scene(); sc.strand("early", [_part("replicated domain: two copies")],
                                label="already replicated")
        sc.strand("late", [_part("unreplicated domain: one copy")],
                                label="not yet replicated")
        sc.note("early", f"one {self.cell_fraction} cell; {self.amplification} preserves relative dosage")
        return sc
