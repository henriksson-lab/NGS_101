"""Hairpin capture of pre-existing DNA ends, as used by END-seq."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Row, Scene, Segment


@dataclass(frozen=True)
class UserHairpin:
    """A single hairpin oligo opened at deoxyuridines by USER enzyme."""
    name: str
    sequence: str
    biotin_dt: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if "U" not in self.sequence:
            raise ValueError(f"{self.name}: USER-openable hairpin needs at least one U")
        if self.sequence != self.sequence.upper():
            raise ValueError(f"{self.name}: store the ordered oligo in uppercase")
        if any(i < 0 or i >= len(self.sequence) or self.sequence[i] != "T"
               for i in self.biotin_dt):
            raise ValueError(f"{self.name}: biotin-dT positions must point to T residues")

    def rows(self) -> list[Row]:
        return [
            Row(chunks=[(f"5'-p {self.sequence} -3'", "r1", False)]),
            Row(chunks=[("folded hairpin: paired Illumina arms joined through a U-bearing loop", None, False)]),
        ]


def captured_break_scene(label: str = "blunted DSB end") -> Scene:
    parts = [Segment(label, "X" * 34, placeholder=True),
             Segment("3-prime dA", "A"),
             Segment("END-seq adaptor 1", "X" * 28, "r1", placeholder=True)]
    sc = Scene.duplex(parts, label="captured genomic end")
    sc.junction("top", "3-prime dA", "END-seq adaptor 1", "hairpin-adaptor ligation")
    sc.labels("top")
    return sc


def user_opening_rows(first: UserHairpin, second: UserHairpin) -> list[Row]:
    """The two closed adapters become PCR-addressable ends only after USER cleavage."""
    opened = [Segment("open P5 arm", "X" * 20, "r1", placeholder=True),
              Segment("break-derived insert", "X" * 34, placeholder=True),
              Segment("open P7 arm", "X" * 20, "r2", placeholder=True)]
    sc = Scene.duplex(opened, label="USER-opened captured molecule")
    sc.junction("top", "open P5 arm", "break-derived insert", "adaptor 1 ligation")
    sc.junction("top", "break-derived insert", "open P7 arm", "adaptor 2 ligation")
    sc.labels("top")
    return [
        Row(chunks=[(f"{first.name}: hairpin -- USER(U) -> open P5 arm", "r1", False)]),
        Row(chunks=[(f"{second.name}: hairpin -- USER(U) -> open P7 arm", "r2", False)]),
        *sc.rows(),
    ]
