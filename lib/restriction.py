"""Restriction digestion, cohesive-end fill-in and ligation for linear DNA.

Unlike :mod:`plasmid`, which finds sites and predicts products on circular maps, this
module keeps the two strand cut coordinates explicit.  It is intended for reaction
schematics: the same ``Digest`` object derives the staggered fragments, their filled
blunt ends, and the sequence created when two ends are ligated.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import complement, revcomp


@dataclass(frozen=True)
class RestrictionEnzyme:
    name: str
    site: str
    top_cut: int       # boundary in site coordinates, 0..len(site)
    bottom_cut: int    # boundary in the bottom strand, drawn left-to-right

    def __post_init__(self) -> None:
        n = len(self.site)
        if not (0 <= self.top_cut <= n and 0 <= self.bottom_cut <= n):
            raise ValueError(f"{self.name}: cut coordinates must lie inside {self.site}")
        if self.top_cut == self.bottom_cut:
            return
        if self.site != revcomp(self.site):
            raise ValueError(f"{self.name}: this symmetric-site model needs a palindrome")

    @property
    def overhang(self) -> str:
        a, b = sorted((self.top_cut, self.bottom_cut))
        return self.site[a:b]

    @property
    def end(self) -> str:
        if self.top_cut == self.bottom_cut:
            return "blunt"
        return "5-prime" if self.top_cut < self.bottom_cut else "3-prime"


MBOI = RestrictionEnzyme("MboI", "GATC", 0, 4)
DPNII = RestrictionEnzyme("DpnII", "GATC", 0, 4)
RSAI = RestrictionEnzyme("RsaI", "GTAC", 2, 2)
ECORV = RestrictionEnzyme("EcoRV", "GATATC", 3, 3)
ECORI = RestrictionEnzyme("EcoRI", "GAATTC", 1, 5)
HINFI = RestrictionEnzyme("HinfI", "GANTC", 1, 4)
BGLII = RestrictionEnzyme("BglII", "AGATCT", 1, 5)
ALUI = RestrictionEnzyme("AluI", "AGCT", 2, 2)


@dataclass(frozen=True)
class FilledDigest:
    digest: "Digest"
    left_new_top: str
    right_new_bottom: str       # 5'->3', not the left-to-right drawn orientation
    biotin_base: str | None = None

    @property
    def junction(self) -> str:
        """Top strand across a ligation between the filled left and right ends."""
        d = self.digest
        return d.enzyme.site[:d.enzyme.bottom_cut] + d.enzyme.site[d.enzyme.top_cut:]

    @property
    def biotin_left_top(self) -> tuple[int, ...]:
        return tuple(i for i, b in enumerate(self.left_new_top) if b == self.biotin_base)

    @property
    def biotin_right_bottom(self) -> tuple[int, ...]:
        return tuple(i for i, b in enumerate(self.right_new_bottom) if b == self.biotin_base)

    @property
    def junction_biotin_top(self) -> tuple[int, ...]:
        """Biotin positions on the ligated top strand, in ``junction`` coordinates."""
        return tuple(self.digest.enzyme.top_cut + i for i in self.biotin_left_top)

    @property
    def junction_biotin_bottom(self) -> tuple[int, ...]:
        """Biotin positions on the drawn bottom strand, in top-strand coordinates."""
        left = self.digest.enzyme.bottom_cut
        n = len(self.right_new_bottom)
        return tuple(left + n - 1 - i for i in self.biotin_right_bottom)


@dataclass(frozen=True)
class Digest:
    """One restriction site with arbitrary flanking top-strand sequence."""
    enzyme: RestrictionEnzyme
    left_flank: str
    right_flank: str

    def __post_init__(self) -> None:
        allowed = set("ACGTNX")
        if not self.left_flank or not self.right_flank:
            raise ValueError("restriction-site drawings need sequence on both sides")
        if set((self.left_flank + self.right_flank).upper()) - allowed:
            raise ValueError("flanks must contain DNA or explicit N/X placeholders")

    @property
    def top(self) -> str:
        return self.left_flank + self.enzyme.site + self.right_flank

    @property
    def bottom_drawn(self) -> str:
        return complement(self.top)

    @property
    def left_top(self) -> str:
        return self.left_flank + self.enzyme.site[:self.enzyme.top_cut]

    @property
    def left_bottom_drawn(self) -> str:
        return complement(self.left_flank + self.enzyme.site[:self.enzyme.bottom_cut])

    @property
    def right_top(self) -> str:
        return self.enzyme.site[self.enzyme.top_cut:] + self.right_flank

    @property
    def right_bottom_drawn(self) -> str:
        return complement(self.enzyme.site[self.enzyme.bottom_cut:] + self.right_flank)

    def fill_in(self, biotin_base: str | None = None) -> FilledDigest:
        """Extend recessed 3' ends to blunt ends.

        ``biotin_base`` names the substituted nucleotide in the polymerase mix (``A``
        for MboI-based in situ Hi-C).  Its positions on each newly made strand are
        derived from the cohesive end.
        """
        if self.enzyme.end != "5-prime":
            raise ValueError("fill_in currently models 5-prime cohesive ends")
        if biotin_base is not None and biotin_base not in "ACGT":
            raise ValueError("biotin_base must be A, C, G or T")
        overhang = self.enzyme.overhang
        return FilledDigest(self, left_new_top=overhang,
                            right_new_bottom=revcomp(overhang),
                            biotin_base=biotin_base)
