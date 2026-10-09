"""Ordered segmented arrays used to increase long-read sequencing throughput."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Row, Scene, Segment


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
        """Draw the complete array as one typed strand with exact junction marks.

        The earlier compact row replaced each insert with prose such as ``[cDNA 1]``.
        That made the SVG renderer see untyped text rather than a molecule and also
        discarded the segments' model-derived hover information.  The linear construct
        already owns the complete ordered sequence, so render that model directly and
        anchor both sides of every covalent segmentation junction by segment name.
        """
        linear = self.linear()
        scene = Scene()
        scene.strand("array", list(linear), label="")
        for i, junction in enumerate(self.junctions):
            left = f"insert {i + 1}: {self.inserts[i].segments[-1].name}"
            right = f"insert {i + 2}: {self.inserts[i + 1].segments[0].name}"
            scene.junction("array", left, junction.name)
            scene.junction("array", junction.name, right)
        return scene.rows()
