"""Reusable states for targeted chromosome-conformation capture.

These constructors encode the distinctions that are easy to lose in prose: a 4C
inverse-PCR template has been cut and ligated twice, whereas Capture-C retains a
linear, adapter-bearing fragment, and a multiway product contains at least two
proximity-ligation junctions in one physical molecule.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment


def region(name: str, width: int = 18) -> Segment:
    return Segment(name, "X" * width, placeholder=True)


@dataclass(frozen=True)
class ThreeCLibrary:
    primary_enzyme: str
    fragments: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.primary_enzyme:
            raise ValueError("3C construction requires a named primary restriction enzyme")
        if len(self.fragments) < 2:
            raise ValueError("a 3C contact requires at least two restriction fragments")

    def scene(self, label: str = "3C proximity-ligation product") -> Scene:
        parts = [region(name) for name in self.fragments]
        sc = Scene(); sc.strand("contact", parts, label=label)
        for left, right in zip(parts, parts[1:]):
            sc.junction("contact", left.name, right.name, "proximity ligation")
        sc.note("contact", f"chromatin digested with {self.primary_enzyme}")
        sc.labels("contact")
        return sc


@dataclass(frozen=True)
class Circular4CTemplate:
    contact: ThreeCLibrary
    secondary_enzyme: str
    viewpoint: str

    def __post_init__(self) -> None:
        if not self.secondary_enzyme:
            raise ValueError("4C circularization requires a named secondary enzyme")
        if self.secondary_enzyme == self.contact.primary_enzyme:
            raise ValueError("4C primary and secondary digests must be distinct")
        if self.viewpoint not in self.contact.fragments:
            raise ValueError("the inverse-PCR viewpoint must occur in the 3C contact")

    def scene(self) -> Scene:
        sc = self.contact.scene("secondary-digested, recircularized 4C template")
        sc.note("contact",
                f"{self.secondary_enzyme} trims the 3C product; dilute ligation closes its ends")
        sc.mark("contact", self.viewpoint, "outward viewpoint primers")
        return sc


def inverse_pcr_product(template: Circular4CTemplate) -> Scene:
    """Open a validated circular template at the outward-facing primer pair."""
    parts = [region("viewpoint-primer end"), region("captured contact DNA"),
             region("viewpoint-primer end 2")]
    sc = Scene(); sc.strand("amplicon", parts, label="viewpoint-centred inverse-PCR amplicon")
    sc.note("amplicon", "amplification proceeds outward from the known viewpoint")
    sc.labels("amplicon")
    return sc


def multiway_contact(primary_enzyme: str, viewpoint: str,
                     reporters: tuple[str, ...]) -> ThreeCLibrary:
    if len(reporters) < 2:
        raise ValueError("multiway capture requires at least two reporter fragments")
    return ThreeCLibrary(primary_enzyme, (viewpoint, *reporters))
