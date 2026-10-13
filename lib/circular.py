"""Single-strand circularisation and rolling-circle amplification.

This is distinct from padlock capture: here an existing linear ssDNA molecule joins its
own 3'-OH to its 5'-phosphate.  The closed template and every RCA repeat are derived from
that one input sequence.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Row, Segment, circle_rows, revcomp, strand_row


@dataclass(frozen=True)
class CircularizationGeometry:
    """Source-scoped topology of a circle-closing reaction.

    A closed product alone cannot reveal whether its ends met directly, on a splint,
    through overlap assembly, or after copying a padlock gap.  Keeping that distinction
    in the reaction model prevents renderers and completeness checks from inferring a
    ligation mechanism merely from the word ``circularize``.
    """

    availability: str
    kind: str
    statement: str

    _KINDS = {
        "single-strand intramolecular ligation",
        "double-strand intramolecular ligation",
        "splinted single-strand ligation",
        "overlap assembly",
        "gap-sealing probe closure",
        "unavailable",
    }

    def __post_init__(self) -> None:
        if self.availability not in {"known", "unavailable"}:
            raise ValueError("circularization geometry must be known or unavailable")
        if self.kind not in self._KINDS:
            raise ValueError(f"unknown circularization geometry {self.kind!r}")
        if not self.statement.strip():
            raise ValueError("circularization geometry needs a source-scoped statement")
        if (self.availability == "unavailable") != (self.kind == "unavailable"):
            raise ValueError("only unavailable evidence may use unavailable geometry")

    @classmethod
    def known(cls, kind: str, *, statement: str) -> "CircularizationGeometry":
        return cls("known", kind, statement)

    @classmethod
    def unavailable(cls, reason: str) -> "CircularizationGeometry":
        return cls("unavailable", "unavailable", reason)


@dataclass(frozen=True)
class ClosedCircle:
    linear: Construct
    closure_left: str
    closure_right: str

    @property
    def sequence(self) -> str:
        return self.linear.top()

    def __len__(self) -> int:
        return len(self.linear)

    def rows(self, closure_label: str = "covalently closed") -> list[Row]:
        """Render the return bond as part of the molecular state, not as prose."""
        return circle_rows(self.linear, closure_label)


def circularize_ssdna(linear: Construct, *, five_prime_phosphate: bool,
                      three_prime_oh: bool = True) -> ClosedCircle:
    """Close a linear ssDNA by joining its last base to its first base."""
    if not len(linear):
        raise ValueError("cannot circularize an empty strand")
    if not five_prime_phosphate:
        raise ValueError("ssDNA circularisation requires a 5'-phosphate")
    if not three_prime_oh:
        raise ValueError("ssDNA circularisation requires a free 3'-OH")
    named = [s.name for s in linear if len(s)]
    if not named or not named[0] or not named[-1]:
        raise ValueError("circle ends need named segments so the closure is unambiguous")
    return ClosedCircle(linear, named[-1], named[0])


@dataclass(frozen=True)
class RCAProduct:
    template: ClosedCircle
    primer: str
    binding_start: int
    unit: str
    copies: int

    @property
    def sequence(self) -> str:
        return self.unit * self.copies

    def __len__(self) -> int:
        return len(self.sequence)

    def construct(self, tag: str = "r3") -> Construct:
        """Return the computed concatemer with one named segment per circle traversal."""
        return Construct([
            Segment(f"circle-complement repeat {number}", self.unit, tag,
                    placeholder=any(base not in "ACGTUN" for base in self.unit))
            for number in range(1, self.copies + 1)
        ], f"{self.copies}-repeat rolling-circle product")

    def rows(self, label: str = "single-stranded rolling-circle concatemer") -> list[Row]:
        """Render the actual computed repeats rather than a prose RCA placeholder."""
        return [strand_row(self.construct(), prefix=f"{label} 5′-", suffix="-3′")]


def rolling_circle(circle: ClosedCircle, primer: str, copies: int = 3) -> RCAProduct:
    """Prime a closed ssDNA and derive tandem complementary copies.

    ``primer`` is written 5'->3'.  It must have one exact site on the circular template.
    One output unit begins with that primer and continues around the rest of the circle;
    therefore its length is guaranteed to equal the template length.
    """
    if copies < 1:
        raise ValueError("RCA must draw at least one copy")
    seq, site = circle.sequence, revcomp(primer)
    if not primer or len(primer) > len(seq):
        raise ValueError("RCA primer must fit on the circular template")
    hay = seq + seq[:len(site) - 1]
    hits = sorted({i % len(seq) for i in range(len(seq)) if hay.startswith(site, i)})
    if len(hits) != 1:
        raise ValueError(f"RCA primer has {len(hits)} sites on {circle.linear.name!r}")
    start = hits[0]
    after = (start + len(site)) % len(seq)
    remainder = "".join(seq[(after + i) % len(seq)]
                        for i in range(len(seq) - len(site)))
    unit = primer + revcomp(remainder)
    if len(unit) != len(seq):
        raise AssertionError("one RCA traversal must copy the circle exactly once")
    return RCAProduct(circle, primer, start, unit, copies)
