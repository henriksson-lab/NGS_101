"""Shared chemistry for DamID-family adenine-marking libraries."""
from __future__ import annotations
from dataclasses import dataclass
from chemdraw import Construct, Segment


@dataclass(frozen=True)
class DamMarkedGATC:
    """A GATC whose adenine methylation is deposited before DNA extraction."""
    target: str
    methylase: str = "E. coli Dam"

    def __post_init__(self) -> None:
        if not self.target:
            raise ValueError("a Dam mark needs a chromatin target")

    def locus(self) -> Construct:
        return Construct([
            Segment("left genomic flank", "X" * 20, placeholder=True),
            Segment("G before m6A", "G", "me"),
            Segment("Dam-deposited m6A", "A", "w1",
                    note=f"methylated near {self.target}"),
            Segment("TC after m6A", "TC", "me"),
            Segment("right genomic flank", "X" * 20, placeholder=True),
        ], name=f"Dam-marked {self.target} locus")

    @property
    def selective_enzyme(self) -> str:
        return "DpnI"


def damid_adapter(top: str, bottom: str) -> Construct:
    """Keep source-transcribed DamID adapter strands coupled as one validated object."""
    if not top or not bottom or "TAATACGACTCACTATA" not in top:
        raise ValueError("DamID adapter top must contain its T7 promoter")
    return Construct([
        Segment("T7-bearing AdRt", top, "t7", bottom=" " * (len(top)-len(bottom)) + bottom),
    ], name="DamID adapter")
