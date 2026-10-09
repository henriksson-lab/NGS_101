"""Ordered segmented arrays used to increase long-read sequencing throughput."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Row, Segment


@dataclass(frozen=True)
class SegmentedArray:
    """A linear array of inserts separated by orientation-specific junctions.

    The constructor owns the alternating insert/junction order.  Callers supply one
    junction for every internal boundary; an incomplete or over-specified array cannot
    be rendered as a complete array accidentally.
    """

    inserts: tuple[Construct, ...]
    junctions: tuple[Segment, ...]
    left_terminal: Segment
    right_terminal: Segment
    name: str = "segmented array"

    def __post_init__(self) -> None:
        if len(self.inserts) < 2:
            raise ValueError("a segmented array needs at least two inserts")
        if len(self.junctions) != len(self.inserts) - 1:
            raise ValueError("a segmented array needs exactly one junction per insert boundary")
        if not len(self.left_terminal) or not len(self.right_terminal):
            raise ValueError("a complete segmented array needs both terminal adapters")
        names = [j.name for j in self.junctions]
        if any(not n for n in names) or len(set(names)) != len(names):
            raise ValueError("segmentation junctions must have unique, non-empty names")

    def core(self) -> Construct:
        """The alternating insert/junction body, before terminal adapters close it."""
        parts = []
        for i, insert in enumerate(self.inserts):
            parts.extend(Segment(f"insert {i + 1}: {s.name}", s.top, s.tag,
                                 s.placeholder, s.inferred, s.bottom, s.note, s.feature)
                         for s in insert)
            if i < len(self.junctions):
                parts.append(self.junctions[i])
        return Construct(parts, name=f"{self.name} core")

    def linear(self) -> Construct:
        return Construct([self.left_terminal, *self.core().segments, self.right_terminal],
                         name=self.name)

    def rows(self) -> list[Row]:
        chunks = [(self.left_terminal.top, self.left_terminal.tag,
                   self.left_terminal.inferred)]
        for i, insert in enumerate(self.inserts):
            label = f"[cDNA {i + 1}]"
            chunks.append((label, "r1", False))
            if i < len(self.junctions):
                chunks.append((" ** ", None, False))
                chunks.append((f"[segment {i + 1}|{i + 2}]",
                               self.junctions[i].tag, self.junctions[i].inferred))
                chunks.append((" ** ", None, False))
        chunks.append((self.right_terminal.top, self.right_terminal.tag,
                       self.right_terminal.inferred))
        return [Row(chunks=chunks)]
