"""Covalently closed dumbbell templates such as PacBio SMRTbell libraries."""
from __future__ import annotations

from dataclasses import dataclass

from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import (Construct, HairpinCapVisual, Row, Segment, StrandVisual,
                      revcomp, strand_row)
from chemdraw import bases_pair
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

    def rows(self) -> list[Row]:
        """Render this validated closed molecule in protocol workflows."""
        return dumbbell_rows(self)


def dumbbell_rows(molecule: Dumbbell) -> list[Row]:
    """Draw the closed topology while retaining the insert's segment colours."""
    gutter = 18
    suffix = " -3'" + " " * 18
    top = strand_row(molecule.insert, "top", prefix=" " * gutter + "5'- ", suffix=suffix)
    bottom = strand_row(molecule.insert, "bottom", prefix=" " * gutter + "3'- ",
                        suffix=" -5'" + " " * 18)
    start = top.visual.start
    end = top.visual.end
    adapter_title = (f"hairpin adapter; {len(molecule.left_adapter)} nt; "
                     f"{molecule.adapter_overhang} overhang pairs with insert "
                     f"{molecule.insert_overhang}")
    top.visual = (
        top.visual,
        HairpinCapVisual(start, "left", "left hairpin", adapter_title),
        HairpinCapVisual(end, "right", "right hairpin", adapter_title),
    )
    # The bottom strand still receives its own directional box; the cap geometry on the
    # top row joins both strand ends one row apart and asserts the molecule has no free end.
    return [top, bottom]


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
