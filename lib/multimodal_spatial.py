"""Reusable molecular views for multimodal and spatial barcode-transfer protocols."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import MolecularFeature, Row, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind


def overload_carriers(graph, template, *, prefix, state, bead_label,
                      inferred=False, source=""):
    """Put intact carriers and barcode beads into deliberately overloaded droplets.

    The rejected droplet population is a mandatory waste terminal.  The returned bead
    remains available to the assay-specific barcode-transfer reaction inside droplets;
    partitioning alone does not pretend that its oligo has become covalent.
    """
    bead = graph.input(
        prefix + "_barcode_bead", bead_label, MaterialKind.OLIGO,
        inferred=inferred,
        source=(source if inferred else None))
    condition = graph.input(
        prefix + "_condition", "overloaded droplet partitioning condition",
        MaterialKind.CONDITION)
    droplets = graph.product(
        prefix + "_contents", "overloaded droplets containing intact carriers and barcode beads",
        MaterialKind.MOLECULAR, state=state, inferred=inferred,
        source=(source if inferred else None))
    waste = graph.product(prefix + "_waste",
                          "droplet material not retained with an intact carrier",
                          MaterialKind.WASTE)
    graph.react(
        prefix + "_partition", "Overload intact carriers into barcode-bead droplets",
        StepKind.PARTITION,
        ((template, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
         (condition, InputRole.CONDITION)), (droplets, waste),
        changes_molecule=False)
    graph.discard(waste, "droplet contents not recovered with an intact carrier")
    reviewed = [bead.id, condition.id, prefix + "_partition"]
    return droplets, bead, reviewed


def recover_and_redistribute_carriers(graph, template, *, prefix,
                                      recovered_state, aliquoted_state):
    """Break droplets, retain intact carriers, then repartition them into wells.

    Cleanup and repartition each expose their rejected material rather than allowing it
    to disappear from the protocol graph.
    """
    recovery_condition = graph.input(
        prefix + "_recovery_condition",
        "emulsion-breaking and intact-carrier recovery condition",
        MaterialKind.CONDITION)
    recovered = graph.product(
        prefix + "_recovered", "recovered intact droplet-barcoded carriers",
        MaterialKind.MOLECULAR, state=recovered_state)
    recovery_waste = graph.product(
        prefix + "_recovery_waste",
        "broken-emulsion material not retained with carriers", MaterialKind.WASTE)
    graph.react(
        prefix + "_recover", "Break droplets and recover intact carriers",
        StepKind.CLEANUP,
        ((template, InputRole.TEMPLATE),
         (recovery_condition, InputRole.CONDITION)),
        (recovered, recovery_waste), changes_molecule=False)
    graph.discard(
        recovery_waste,
        "emulsion and bead material removed during intact-carrier recovery")

    aliquot_condition = graph.input(
        prefix + "_aliquot_condition", "aliquot redistribution condition",
        MaterialKind.CONDITION)
    aliquoted = graph.product(
        prefix + "_aliquoted",
        "droplet-barcoded carriers redistributed to indexed aliquots",
        MaterialKind.MOLECULAR, state=aliquoted_state)
    aliquot_waste = graph.product(
        prefix + "_aliquot_waste",
        "carrier material outside the retained aliquot representation",
        MaterialKind.WASTE)
    graph.react(
        prefix + "_redistribute",
        "Mix recovered carriers and redistribute them into indexed aliquots",
        StepKind.PARTITION,
        ((recovered, InputRole.TEMPLATE),
         (aliquot_condition, InputRole.CONDITION)),
        (aliquoted, aliquot_waste), changes_molecule=False)
    graph.discard(
        aliquot_waste,
        "not represented among retained aliquot-indexed carrier molecules")
    reviewed = [
        recovery_condition.id, prefix + "_recover", aliquot_condition.id,
        prefix + "_redistribute",
    ]
    return aliquoted, reviewed


def finish_pcr_branch(graph, prefix, template, final, primers, reviewed, *, label,
                      read_lengths=None, read_uses=None, inferred_primers=False,
                      limitation="The source defines the finished library but not one universal PCR-primer pair.",
                      geometry=None):
    """PCR-complete and sequence one branch without duplicating endpoint logic.

    The caller remains responsible for every upstream molecular operation.  This helper
    centralizes the invariant common to multimodal protocols: each physical branch gets
    its own primer inputs, typed PCR step, library terminal and primer-derived readout.
    """
    source = graph.authority + (" (endpoint primer sequence not declared in this model)" if inferred_primers else "")
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO, inferred=inferred_primers, source=source)
    reverse = graph.input(prefix + "_reverse", label + " reverse PCR primer",
                          MaterialKind.OLIGO, inferred=inferred_primers, source=source)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final,
                            inferred=inferred_primers, source=source)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " branch",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=(geometry or PCRGeometry.unavailable(
                    final, reason=limitation)))
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=tuple(p.role for p in primers),
        read_lengths=read_lengths or {}, read_uses=read_uses or {}))
    reviewed.extend((forward.id, reverse.id, polymerase.id,
                     prefix + "_pcr", library.id))
    return library


@dataclass(frozen=True)
class HybridizationLigationRound:
    """One split-pool barcode oligo pre-annealed to a two-arm splint.

    The splint must be exactly the reverse complement of the prior molecule's
    exposed 5′ arm followed by the reverse complement of the barcode oligo's
    ligating 3′ arm.  This makes both overhangs and the covalent product order
    consequences of the ordered oligos rather than separately drawn claims.
    """

    name: str
    previous_5_arm: str
    barcode_5_arm: str
    barcode_3_arm: str
    linker: str
    feature: MolecularFeature
    phosphorylated_5: bool = False
    barcode_length: int = 8

    def __post_init__(self) -> None:
        expected = revcomp(self.previous_5_arm) + revcomp(self.barcode_3_arm)
        if self.linker != expected:
            raise ValueError(
                f"{self.name}: linker does not bridge the prior molecule and barcode oligo"
            )
        if self.barcode_length < 1:
            raise ValueError(f"{self.name}: barcode length must be positive")

    def product_segments(self) -> list[Segment]:
        """The barcode strand in 5′→3′ product order."""
        return [
            Segment(f"{self.name} outgoing arm", self.barcode_5_arm),
            Segment(f"{self.name} barcode", "B" * self.barcode_length, "cbc",
                    placeholder=True, feature=self.feature),
            Segment(f"{self.name} ligating arm", self.barcode_3_arm),
        ]

    @property
    def blocking_oligo(self) -> str:
        """Sequence that sequesters an unbound barcode/linker assembly."""
        return revcomp(self.linker)

    def adapter_scene(self) -> Scene:
        """The complete splinted nick, including both molecules being joined.

        The preceding molecule ends in ``previous_5_arm`` and supplies the
        obligatory 5′ phosphate.  The new barcode oligo supplies the 3′ OH at
        the end of ``barcode_3_arm``.  Both anneal next to one another on the
        linker; omitting the preceding strand would turn the actual ligation
        partner into an unexplained label on the splint.
        """
        linker = [
            Segment(f"{self.name} prior-molecule splint", revcomp(self.previous_5_arm)),
            Segment(f"{self.name} barcode splint", revcomp(self.barcode_3_arm)),
        ]
        barcode = self.product_segments()
        sc = Scene()
        sc.strand(f"{self.name} linker", linker, label=f"{self.name} linker")
        sc.anneal(
            f"{self.name} preceding molecule",
            [Segment(f"{self.name} preceding 5-prime arm", self.previous_5_arm)],
            to=f"{self.name} linker",
            pair=(f"{self.name} preceding 5-prime arm",
                  f"{self.name} prior-molecule splint"),
            label=f"{self.name} preceding molecule", mod5="/5Phos/",
        )
        sc.anneal(
            f"{self.name} barcode oligo", barcode,
            to=f"{self.name} linker",
            pair=(f"{self.name} ligating arm", f"{self.name} barcode splint"),
            label=f"{self.name} barcode oligo",
            mod5="/5Phos/" if self.phosphorylated_5 else "",
        )
        sc.ligation_end(
            f"{self.name} preceding molecule",
            f"{self.name} preceding 5-prime arm",
            kind="splinted", polarity="5-prime", owner="substrate",
            text="5′-phosphate at the ligation nick",
        )
        sc.ligation_end(
            f"{self.name} barcode oligo", f"{self.name} ligating arm",
            kind="splinted", polarity="3-prime", owner="adapter",
            text="3′-OH at the ligation nick",
        )
        sc.labels(f"{self.name} linker")
        sc.labels(f"{self.name} preceding molecule")
        sc.labels(f"{self.name} barcode oligo")
        return sc


def barcode_rounds(rounds: int, *, label: str = "cell barcode") -> list[Row]:
    """Draw split-pool barcode ligations, with every covalent junction explicit."""
    if rounds < 1:
        raise ValueError("at least one barcode round is required")
    text = "molecule"
    for n in range(1, rounds + 1):
        text += f" ** [{label} {n}]"
    return [Row(chunks=[(text, "cbc", False)]),
            Row(chunks=[("** = ligation at the named barcode boundary", "me", False)])]


def modality_split(*names: str) -> list[Row]:
    """Show one tagged cell yielding separately amplified sequencing libraries."""
    if len(names) < 2:
        raise ValueError("a modality split needs at least two libraries")
    rows = [Row(chunks=[("one cell identity", "cbc", False)])]
    for i, name in enumerate(names):
        rows.append(Row(chunks=[(("├─ " if i < len(names) - 1 else "└─ ") + name,
                                 None, False)]))
    return rows


def adjacent_probe_scene() -> Scene:
    """Adjacent probe halves on fixed RNA; Scene enforces the duplex geometry."""
    target = [Segment("target RNA left", "X" * 18, placeholder=True),
              Segment("target RNA right", "X" * 18, placeholder=True)]
    left = [Segment("left probe barcode", "L" * 12, "cbc", placeholder=True,
                    feature=feature("probe_barcode_left", "feature_barcode", "unknown",
                                    group="probe_barcode", part="left")),
            Segment("left probe target arm", "X" * 18, placeholder=True)]
    right = [Segment("right probe target arm", "X" * 18, placeholder=True),
             Segment("right probe barcode", "R" * 12, "cbc", placeholder=True,
                     feature=feature("probe_barcode_right", "feature_barcode", "unknown",
                                     group="probe_barcode", part="right"))]
    sc = Scene(); sc.strand("RNA", target, label="fixed target RNA")
    sc.anneal("left probe", left, to="RNA", pair=("left probe target arm", "target RNA left"),
              label="left probe", unpaired=("left probe barcode",))
    sc.anneal("right probe", right, to="RNA", pair=("right probe target arm", "target RNA right"),
              label="right probe", unpaired=("right probe barcode",))
    sc.arrow("left probe", "adjacent probe-pair ligation")
    sc.labels("left probe"); sc.labels("right probe")
    return sc


def surface_capture_rows(surface: str, barcode: str, molecule: str) -> list[Row]:
    return [
        Row(chunks=[(f"{surface}:  ●—[{barcode}]—[UMI]—poly(dT)", "cbc", False)]),
        Row(chunks=[(f"                                      || {molecule} poly(A)", None, False)]),
        Row(chunks=[("                                      ← reverse transcription", "me", False)]),
    ]
