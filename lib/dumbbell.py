"""Covalently closed dumbbell templates such as PacBio SMRTbell libraries."""
from __future__ import annotations

from dataclasses import dataclass

from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import (Construct, HairpinCapVisual, Row, Segment, StrandVisual,
                      YAdapterVisual, revcomp, strand_row)
from chemdraw import bases_pair
from hairpin import OpenedHairpin
from protocol_graph import SequencingPlan


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
    left_stem_length: int = 0
    right_stem_length: int = 0

    def __post_init__(self) -> None:
        if not len(self.insert):
            raise ValueError("a dumbbell needs a duplex insert")
        if not self.left_adapter or not self.right_adapter:
            raise ValueError("both insert ends need a hairpin adapter")
        if bool(self.insert_overhang) != bool(self.adapter_overhang):
            raise ValueError("insert and hairpin must either both declare an overhang or neither")
        if self.insert_overhang and not bases_pair(self.insert_overhang,
                                                   self.adapter_overhang):
            raise ValueError("hairpin and insert overhangs are not complementary")
        for side, adapter, stem in (("left", self.left_adapter, self.left_stem_length),
                                    ("right", self.right_adapter, self.right_stem_length)):
            if stem < 0 or 2 * stem + len(self.adapter_overhang) > len(adapter):
                raise ValueError(f"{side} hairpin has invalid stem/overhang lengths")
            if stem:
                end = len(adapter) - len(self.adapter_overhang)
                if adapter[:stem] != revcomp(adapter[end - stem:end]):
                    raise ValueError(f"{side} hairpin stem does not self-pair")

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

    def rows(self) -> list[Row]:
        """Render this validated closed molecule in protocol workflows."""
        return dumbbell_rows(self)


def dumbbell_rows(molecule: Dumbbell) -> list[Row]:
    """Draw the closed topology while retaining the insert's segment colours."""
    left_extent = max(28, len(molecule.left_adapter) // 2
                      + molecule.left_stem_length + 4)
    right_extent = max(28, len(molecule.right_adapter) // 2
                       + molecule.right_stem_length + 4)
    gutter = left_extent
    suffix = " -3'" + " " * right_extent
    top = strand_row(molecule.insert, "top", prefix=" " * gutter + "5'- ", suffix=suffix)
    bottom = strand_row(molecule.insert, "bottom", prefix=" " * gutter + "3'- ",
                        suffix=" -5'" + " " * right_extent)
    start = top.visual.start
    end = top.visual.end
    def adapter_title(side: str, sequence: str) -> str:
        title = Segment(f"{side} hairpin adapter", sequence).hover_text()
        if molecule.adapter_overhang:
            return (title + f"; {molecule.adapter_overhang} overhang pairs with insert "
                    + molecule.insert_overhang)
        return title + "; direct repaired junction to insert"
    top.visual = (
        top.visual,
        HairpinCapVisual(start, "left", "left hairpin", molecule.left_adapter,
                         molecule.left_stem_length, len(molecule.adapter_overhang),
                         adapter_title("left", molecule.left_adapter)),
        HairpinCapVisual(end, "right", "right hairpin", molecule.right_adapter,
                         molecule.right_stem_length, len(molecule.adapter_overhang),
                         adapter_title("right", molecule.right_adapter)),
    )
    # The bottom strand still receives its own directional box; the cap geometry on the
    # top row joins both strand ends one row apart and asserts the molecule has no free end.
    # Routed loops and their base labels extend above/below the two duplex rows. Blank
    # rows reserve that vertical space in both standalone panels and reaction workflows.
    return [Row(), Row(), top, bottom, Row(), Row()]


@dataclass(frozen=True)
class OpenedDumbbell:
    """A two-hairpin library after both loops have been cut into Y adapters."""

    closed: Dumbbell
    left: OpenedHairpin
    right: OpenedHairpin
    name: str = "USER-opened Y-adapter library"

    def __post_init__(self) -> None:
        for side, adapter, opened in (
                ("left", self.closed.left_adapter, self.left),
                ("right", self.closed.right_adapter, self.right)):
            rebuilt = opened.left_arm + opened.removed_base + opened.right_arm
            if adapter != rebuilt:
                raise ValueError(f"{side} opened arms do not reconstruct the closed hairpin")

    @staticmethod
    def _title(name: str, sequence: str, end: str) -> str:
        return Segment(name, sequence).hover_text() + f"; new {end} end"

    def rows(self) -> list[Row]:
        longest_left = max(len(self.left.left_arm), len(self.left.right_arm))
        longest_right = max(len(self.right.left_arm), len(self.right.right_arm))
        gutter = longest_left + 12
        tail = " " * (longest_right + 12)
        top = strand_row(self.closed.insert, "top", prefix=" " * gutter + "5'- ",
                         suffix=" -3'" + tail)
        bottom = strand_row(self.closed.insert, "bottom", prefix=" " * gutter + "3'- ",
                            suffix=" -5'" + tail)
        top.visual = (
            top.visual,
            YAdapterVisual(
                top.visual.start, "left", "Read 1 arm",
                self.left.right_arm,
                self._title("left Read 1 arm", self.left.right_arm,
                            self.left.right_new_end),
                "Read 2 arm",
                self.left.left_arm,
                self._title("left Read 2 arm", self.left.left_arm,
                            self.left.left_new_end),
                self.closed.left_stem_length,
            ),
            YAdapterVisual(
                top.visual.end, "right", "Read 2 arm",
                self.right.left_arm,
                self._title("right Read 2 arm", self.right.left_arm,
                            self.right.left_new_end),
                "Read 1 arm",
                self.right.right_arm,
                self._title("right Read 1 arm", self.right.right_arm,
                            self.right.right_new_end),
                self.closed.right_stem_length,
            ),
        )
        return [Row(), Row(), top, bottom, Row(), Row()]


def smrtbell_circuit(molecule: Dumbbell) -> Construct:
    """Build one checked polymerase circuit from a closed dumbbell topology.

    The caller supplies the insert and both hairpins once through :class:`Dumbbell`;
    sequencing plans and protocol graphs therefore cannot silently disagree about the
    order or orientation traversed by the polymerase.
    """
    circuit = Construct([
        Segment("insert strand 1", molecule.insert.top(), placeholder=True,
                length_bp=len(molecule.insert)),
        Segment("right hairpin / primer site", molecule.right_adapter, "r1",
                placeholder=True, length_bp=len(molecule.right_adapter)),
        Segment("insert strand 2 reverse complement", revcomp(molecule.insert.top()),
                placeholder=True, length_bp=len(molecule.insert)),
        Segment("left hairpin / primer site", molecule.left_adapter, "r1",
                placeholder=True, length_bp=len(molecule.left_adapter)),
    ], name=f"one {molecule.name} polymerase circuit")
    if circuit.top() != molecule.pass_sequence:
        raise AssertionError("polymerase circuit must follow the closed dumbbell")
    return circuit


def smrtbell_hairpin_geometry(molecule: Dumbbell, *, statement: str
                              ) -> AdapterLigationGeometry:
    """Describe the shared cohesive end enforced by a :class:`Dumbbell`.

    Hairpin bodies may remain proprietary, but the constructor has already checked that
    the insert and adapter overhangs pair.  Rendering and graph validation can therefore
    share one source of truth for the ligating terminal base.
    """
    return reported_adapter_end(
        f"{molecule.name} hairpin", kind="cohesive", polarity="3-prime",
        terminal_sequence=molecule.adapter_overhang, statement=statement)


def smrtbell_sequencing_plan(molecule: Dumbbell, circuit: Construct, *,
                             limitation: str) -> SequencingPlan:
    """Declare a primer-sequence-limited PacBio endpoint with checked topology."""
    if circuit.top() != molecule.pass_sequence:
        raise ValueError("a PacBio sequencing circuit must follow the closed dumbbell")
    if not limitation.strip():
        raise ValueError("an undisclosed SMRTbell primer needs a source limitation")
    return SequencingPlan(
        library=circuit, primers=(), acquired_roles=("Polymerase read",),
        read_uses={"Polymerase read": "include HiFi consensus from repeated circuits"},
        limitation=limitation)
