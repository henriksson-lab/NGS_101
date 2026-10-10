"""Reusable molecular models for contact assays that do not use ordinary Hi-C ligation."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Segment


@dataclass(frozen=True)
class BivalentTransposaseLinker:
    """A two-ended Tn5 linker constrained to bridge two chromatin targets.

    The constructor owns the defining physical invariant: two complete mosaic ends must
    flank a non-empty spacer.  Protocol modules may then decorate the spacer (for example
    with biotin) without reimplementing the topology.
    """

    mosaic_end: str
    spacer: str
    name: str = "bivalent ME linker"

    def __post_init__(self) -> None:
        if len(self.mosaic_end) != 19:
            raise ValueError("a Tn5 mosaic end must contain 19 nucleotides")
        if not self.spacer:
            raise ValueError("a bivalent transposase linker needs a spacer")

    def linker(self) -> Construct:
        return Construct([
            Segment("ME A", self.mosaic_end, "me"),
            Segment("linker spacer", self.spacer, "linker"),
            Segment("ME B", self.mosaic_end, "me"),
        ], name=self.name)

    def bridged_product(self) -> Construct:
        return Construct([
            Segment("contact locus A", "X" * 24, placeholder=True),
            *list(self.linker()),
            Segment("contact locus B", "X" * 24, placeholder=True),
        ], name=f"{self.name} trans-contact product")
