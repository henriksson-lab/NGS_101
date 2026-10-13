"""Reusable structural model for adjacent-probe branched-DNA RNA detection assays."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Row, Scene, Segment, revcomp


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


@dataclass(frozen=True)
class SequenceExample:
    """A disclosed implementation without implying identity to a current product."""
    title: str
    note: str
    oligos: tuple[tuple[str, tuple[Segment, ...]], ...]
    scene: Scene
    caption: str


def paired_tail_docking_scene(left_tail: str, right_tail: str,
                              preamplifier_site: str) -> Scene:
    """Draw and validate a preamplifier across two adjacent probe-tail halves.

    The two top segments are separate probe molecules conceptually; joining them in one
    Scene strand supplies a common coordinate system for the 28-base docking duplex.
    The footer makes that representational boundary explicit.
    """
    tails = left_tail + right_tail
    expected = revcomp(tails)
    if preamplifier_site.upper() != expected:
        raise ValueError("preamplifier site is not the reverse complement of both Z tails")
    sc = Scene()
    sc.strand("paired Z tails", [
        Segment("left-probe Z tail", left_tail, "cbc"),
        Segment("right-probe Z tail", right_tail, "cbc"),
    ], label="two adjacent Z-probe tails")
    sc.anneal("preamplifier", [
        Segment("preamplifier docking site", preamplifier_site, "r1"),
    ], to="paired Z tails",
              pair=("preamplifier docking site", "left-probe Z tail"),
              label="preamplifier", above=True)
    sc.mark("paired Z tails", "left-probe Z tail", "tail supplied by probe 1")
    sc.mark("paired Z tails", "right-probe Z tail", "tail supplied by probe 2")
    sc.footer("the two tail segments belong to separate target probes",
              "paired Z tails")
    return sc


def adjacent_target_example_scene(left_target: str, right_target: str) -> Scene:
    """Draw disclosed adjacent target regions with the implied probe-binding arms."""
    sc = Scene()
    sc.strand("target RNA", [
        Segment("left disclosed target region", left_target, "me"),
        Segment("right disclosed target region", right_target, "me"),
    ], label="published target-RNA regions")
    sc.anneal("left probe", [
        Segment("left universal tail", "[undisclosed tail]", "cbc", placeholder=True),
        Segment("left target-binding arm", revcomp(left_target), "r1"),
    ], to="target RNA", pair=("left target-binding arm", "left disclosed target region"),
              label="implied left probe arm", above=True,
              unpaired=("left universal tail",))
    sc.anneal("right probe", [
        Segment("right target-binding arm", revcomp(right_target), "r2"),
        Segment("right universal tail", "[undisclosed tail]", "cbc", placeholder=True),
    ], to="target RNA", pair=("right target-binding arm", "right disclosed target region"),
              label="implied right probe arm", above=True,
              unpaired=("right universal tail",))
    sc.footer("binding arms are derived complements; universal tails were not disclosed",
              "target RNA")
    return sc


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
