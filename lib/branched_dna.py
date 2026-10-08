"""Reusable structural model for adjacent-probe branched-DNA RNA detection assays."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Row, Scene, Segment


@dataclass(frozen=True)
class BranchModel:
    """Published multiplicities for one target RNA and one amplification tree."""
    target_pairs_min: int
    target_pairs_max: int
    amplifiers_per_preamp: int
    labels_per_amplifier: int

    def __post_init__(self) -> None:
        values = (self.target_pairs_min, self.target_pairs_max,
                  self.amplifiers_per_preamp, self.labels_per_amplifier)
        if any(n <= 0 for n in values) or self.target_pairs_min > self.target_pairs_max:
            raise ValueError(f"invalid branched-DNA multiplicities: {values}")

    @property
    def labels_per_tree(self) -> int:
        return self.amplifiers_per_preamp * self.labels_per_amplifier

    @property
    def labels_per_target(self) -> tuple[int, int]:
        return (self.target_pairs_min * self.labels_per_tree,
                self.target_pairs_max * self.labels_per_tree)


def target_scene() -> Scene:
    """One adjacent probe pair, placed by the shared strand-pairing machinery.

    X is an undisclosed base, not a claimed sequence.  The widths are schematic; the
    protocol-specific caption reports any published length range.
    """
    target = [
        Segment("left target site", "X" * 12, placeholder=True),
        Segment("right target site", "X" * 12, placeholder=True),
    ]
    left_probe = [
        Segment("left Z tail", "X" * 6, "cbc", placeholder=True),
        Segment("left probe", "X" * 12, "me", placeholder=True),
    ]
    right_probe = [
        Segment("right probe", "X" * 12, "me", placeholder=True),
        Segment("right Z tail", "X" * 6, "cbc", placeholder=True),
    ]
    sc = Scene()
    sc.strand("target RNA", target, label="target RNA")
    sc.anneal("left Z", left_probe, to="target RNA",
              pair=("left probe", "left target site"),
              label="left target probe", above=True,
              unpaired=("left Z tail",))
    sc.anneal("right Z", right_probe, to="target RNA",
              pair=("right probe", "right target site"),
              label="right target probe", above=True,
              unpaired=("right Z tail",))
    sc.mark("left Z", "left Z tail", "half of preamplifier landing site")
    sc.mark("right Z", "right Z tail", "half of preamplifier landing site")
    sc.footer("schematic widths; X = undisclosed base", "target RNA")
    return sc


def tree_rows(model: BranchModel) -> list[Row]:
    """One representative signal tree, annotated with model-derived multiplicities."""
    a, l = model.amplifiers_per_preamp, model.labels_per_amplifier
    return [
        Row(chunks=[(f"label probes       {l} per amplifier", "umi", True)]),
        Row(chunks=[("                    | | | | |", None, False)]),
        Row(chunks=[(f"amplifiers         {a} branches per preamplifier", "r2", True)]),
        Row(chunks=[("                    \\ | | /", None, False)]),
        Row(chunks=[("preamplifier   ===== signal-amplification trunk =====", "r1", True)]),
        Row(chunks=[("                          ||", None, False)]),
        Row(chunks=[("paired tails       [Z-left][Z-right]", "cbc", True)]),
        Row(chunks=[("                          ||", None, False)]),
        Row(chunks=[("target RNA     5'-...adjacent target region...-3'", None, False)]),
    ]
