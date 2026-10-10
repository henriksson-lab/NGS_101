"""End-state gate for nuclease-directed selective adapter ligation.

Methods such as nCATS first remove phosphates from every pre-existing DNA end, then
create fresh nuclease ends that retain the 5'-phosphate required by ligase.  Keeping that
selectivity in a constructor prevents a protocol page from accidentally depicting old and
new ends as equally ligatable.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct


@dataclass(frozen=True)
class SelectiveLigationProduct:
    substrate: Construct
    selected_fragment: Construct
    old_ends_phosphorylated: bool
    fresh_ends_phosphorylated: bool

    def __post_init__(self) -> None:
        if not len(self.substrate) or not len(self.selected_fragment):
            raise ValueError("selective ligation needs non-empty substrate and selected DNA")
        if self.old_ends_phosphorylated:
            raise ValueError("pre-existing ends must be dephosphorylated before selection")
        if not self.fresh_ends_phosphorylated:
            raise ValueError("fresh nuclease ends must retain ligatable 5-prime phosphate")


def dephosphorylate_then_cut(substrate: Construct,
                             selected_fragment: Construct) -> SelectiveLigationProduct:
    """Build the only valid phosphate state for fresh-cut selective ligation."""
    return SelectiveLigationProduct(substrate, selected_fragment, False, True)
