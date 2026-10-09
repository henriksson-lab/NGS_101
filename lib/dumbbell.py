"""Covalently closed dumbbell templates such as PacBio SMRTbell libraries."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Row, revcomp, strand_row
from chemdraw import bases_pair


@dataclass(frozen=True)
class Dumbbell:
    """A duplex insert whose two strand ends are joined by hairpin adapters.

    ``left_adapter`` and ``right_adapter`` are each written in the direction encountered
    while traversing the closed strand: top insert -> right hairpin -> reverse-complement
    insert -> left hairpin.  The resulting molecule has no free DNA ends.
    """
    insert: Construct
    left_adapter: str
    right_adapter: str
    name: str = "dumbbell library"
    insert_overhang: str = "A"
    adapter_overhang: str = "T"

    def __post_init__(self) -> None:
        if not len(self.insert):
            raise ValueError("a dumbbell needs a duplex insert")
        if not self.left_adapter or not self.right_adapter:
            raise ValueError("both insert ends need a hairpin adapter")
        if not bases_pair(self.insert_overhang, self.adapter_overhang):
            raise ValueError("hairpin and insert overhangs are not complementary")

    @property
    def circular_sequence(self) -> str:
        return (self.insert.top() + self.right_adapter
                + revcomp(self.insert.top()) + self.left_adapter)

    @property
    def pass_sequence(self) -> str:
        """One polymerase circuit: forward insert then reverse-complement insert."""
        return self.circular_sequence

    def polymerase_product(self, circuits: int = 2) -> str:
        if circuits < 1:
            raise ValueError("polymerase must make at least one circuit")
        return self.pass_sequence * circuits


def dumbbell_rows(molecule: Dumbbell) -> list[Row]:
    """Draw the closed topology while retaining the insert's segment colours."""
    left, right = "left hairpin", "right hairpin"
    gutter = len(left) + 2
    between = len(molecule.insert) + 14
    # Keep the hairpin/junction punctuation outside the molecular span, but construct
    # both insert strands through the shared typed-row helper.  This gives SVG renderers
    # direction, segment hover text and computed Tm without duplicating that logic here.
    top = strand_row(molecule.insert, "top",
                     prefix=" " * gutter + "|** 5'- ", suffix=" -3' **|")
    bottom = strand_row(molecule.insert, "bottom",
                        prefix=" " * gutter + "|** 3'- ", suffix=" -5' **|")
    cap = Row(chunks=[(left + "  ." + " " * between + ".  " + right,
                       None, False)])
    close = Row(chunks=[(" " * gutter + "'" + "-" * between + "'",
                         None, False)])
    return [cap, top, bottom, close]
