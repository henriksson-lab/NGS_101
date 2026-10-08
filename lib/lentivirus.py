"""Structural transformations of a lentiviral transfer genome.

The two LTRs of the packaged RNA are deliberately asymmetric: the 5' end supplies
R--U5 and the 3' end supplies U3--R. Reverse transcription rebuilds both proviral
LTRs as U3--R--U5. Consequently, cargo placed in the 3' U3 is copied into both
integrated LTRs.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

from chemdraw import Construct, Segment


def _named(prefix: str, segments: Sequence[Segment]) -> list[Segment]:
    """Copy segments while making every name unique in the resulting construct."""
    return [replace(s, name=f"{prefix} {s.name}" if s.name else prefix) for s in segments]


@dataclass(frozen=True)
class TransferGenome:
    """The four parts needed to derive packaged RNA and an integrated provirus.

    ``u3`` is supplied once, from the 3' LTR. :meth:`provirus` necessarily uses that
    same tuple at both ends; callers cannot provide two independently drifting copies.
    """

    u3: tuple[Segment, ...]
    r: tuple[Segment, ...]
    u5: tuple[Segment, ...]
    internal: tuple[Segment, ...]
    name: str = "lentiviral transfer genome"

    def __post_init__(self) -> None:
        if not self.u3 or not self.r or not self.u5:
            raise ValueError("a transfer genome needs non-empty U3, R and U5 regions")

    def viral_rna(self) -> Construct:
        """Packaged 5'-R--U5--internal--U3--R-3' RNA, drawn as DNA bases."""
        return Construct([
            *_named("5' LTR", self.r), *_named("5' LTR", self.u5),
            *_named("internal", self.internal), *_named("3' LTR", self.u3),
            *_named("3' LTR", self.r),
        ], name=f"{self.name} packaged RNA")

    def provirus(self) -> Construct:
        """Integrated U3--R--U5 ... U3--R--U5 provirus derived from this genome."""
        return Construct([
            *_named("5' LTR copied", self.u3), *_named("5' LTR", self.r),
            *_named("5' LTR", self.u5), *_named("internal", self.internal),
            *_named("3' LTR", self.u3), *_named("3' LTR", self.r),
            *_named("3' LTR copied", self.u5),
        ], name=f"{self.name} integrated provirus")
