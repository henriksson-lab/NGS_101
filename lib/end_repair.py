"""Sequence-derived repair of heterogeneous 5-prime DNA overhangs."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import is_real_dna, revcomp


@dataclass(frozen=True)
class FilledOverhang:
    """A 5′ overhang and the strand polymerase must synthesize to make it blunt."""
    overhang_5p: str
    fill_5p: str
    substituted_bases: frozenset[str]

    @property
    def substituted_positions(self) -> tuple[int, ...]:
        """Positions in the newly synthesized 5′→3′ strand carrying analogues."""
        return tuple(i for i, base in enumerate(self.fill_5p)
                     if base in self.substituted_bases)


def fill_five_prime_overhang(overhang_5p: str,
                             substituted_bases=()) -> FilledOverhang:
    """Derive polymerase fill-in for a heterogeneous 5′ overhang.

    ``substituted_bases`` names bases supplied as labelled analogues, such as A and C in
    Micro-C.  The sequence is derived here so diagrams cannot hand-enter an incompatible
    fill strand.
    """
    overhang = overhang_5p.upper()
    labels = frozenset(b.upper() for b in substituted_bases)
    if not is_real_dna(overhang) or set(overhang) - set("ACGT"):
        raise ValueError("5-prime overhang must be non-empty unambiguous DNA")
    if labels - set("ACGT"):
        raise ValueError("substituted bases must be A, C, G or T")
    return FilledOverhang(overhang, revcomp(overhang), labels)
