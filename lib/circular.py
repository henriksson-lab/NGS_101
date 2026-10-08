"""Single-strand circularisation and rolling-circle amplification.

This is distinct from padlock capture: here an existing linear ssDNA molecule joins its
own 3'-OH to its 5'-phosphate.  The closed template and every RCA repeat are derived from
that one input sequence.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, revcomp


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
