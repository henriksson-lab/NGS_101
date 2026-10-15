"""Typed, provenance-preserving process graph for library-preparation protocols.

``Workflow`` is useful for drawing one linear reaction path.  A protocol is broader:
several materials may enter a reaction, a state may split into parallel branches, and a
reaction may yield a retained library plus material that is discarded.  ``ProtocolGraph``
models that process as a directed acyclic *hypergraph*: a step consumes one or more
material nodes and creates one or more new material nodes.

The graph deliberately separates three questions:

* construction makes duplicate producers, backward references and cycles impossible;
* typed step validators establish the local shape of PCR, ligation, digestion, etc.;
* ``publication_problems`` checks global provenance, terminal accounting, sequencing
  readiness and coverage of the authoritative source.

It does not claim that an omitted source step can be discovered computationally.  Source
items are therefore explicit inputs to the completeness contract, alongside the molecular
graph.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import re
from typing import Any, Iterable, Mapping


class MaterialKind(str, Enum):
    BIOLOGICAL = "biological sample"
    MOLECULAR = "molecular material"
    OLIGO = "oligo"
    ADAPTER = "adapter"
    LIBRARY = "sequencing library"
    CATALYST = "enzyme or catalyst"
    SELECTION = "selection reagent"
    CONDITION = "condition or buffer"
    WASTE = "discarded material"
    QC = "QC material"


class InputRole(str, Enum):
    TEMPLATE = "template"
    STRUCTURAL = "structural reagent"
    CATALYST = "catalyst"
    SELECTION = "selection reagent"
    ANALYTE = "measured analyte"
    CONDITION = "condition"
    SPIKE_IN = "spike-in"
    SELF_PRIMER = "intramolecular primer"


class StepKind(str, Enum):
    GENERIC = "generic transformation"
    MAGIC = "magic happens"
    INTEGRATION = "genomic integration"
    GENOME_EDITING = "genome editing"
    CROSS_LINKING = "cross-linking"
    PERMEABILIZATION = "permeabilization"
    SOLUBILIZATION = "solubilization"
    QUENCHING = "reaction quenching"
    THERMAL_TREATMENT = "thermal treatment"
    FRAGMENTATION = "fragmentation"
    EXTRACTION = "extraction"
    DEPROTEINIZATION = "deproteinization"
    END_REPAIR = "end repair"
    REVERSE_TRANSCRIPTION = "reverse transcription"
    TEMPLATE_SWITCHING = "template switching"
    TRANSCRIPTION = "transcription"
    RUN_ON = "nuclear run-on"
    PCR = "PCR"
    LIGATION = "ligation"
    PROXIMITY_LIGATION = "proximity ligation"
    DIGESTION = "digestion"
    TAGMENTATION = "tagmentation"
    ADAPTER_ATTACHMENT = "non-ligase adapter attachment"
    HYBRIDIZATION = "hybridization"
    BINDING = "binding"
    FOLDING = "intramolecular folding"
    STRAND_SEPARATION = "strand separation"
    EXTENSION = "extension or fill-in"
    NICK_SEALING = "nick sealing"
    TAILING = "tailing"
    PHOSPHORYLATION = "phosphorylation"
    DEPHOSPHORYLATION = "dephosphorylation"
    CHEMICAL_CONVERSION = "chemical conversion"
    CIRCULARIZATION = "circularization"
    SELECTION = "selection"
    SAMPLING = "sampling or bottleneck"
    AFFINITY_CAPTURE = "affinity capture"
    CLEANUP = "cleanup"
    SPLIT = "split"
    POOL = "pool"
    PARTITION = "partition"


class MolecularRelation(str, Enum):
    """A narrow molecular claim to check, independent of the displayed operation."""
    CIRCLE_CLOSURE = "circle closure preserves strand order"
    RESTRICTION_FILL = "restriction end fill follows cut coordinates"
    ROLLING_CIRCLE_COPY = "rolling-circle repeats copy the circle"
    CONTACT_JUNCTION = "contact product contains the derived junction"
    HAIRPIN_STEM = "intramolecular stem can pair"
    STRAND_ORDER = "covalent strand order is preserved"
    END_COMBINATIONS = "end combinations have the declared primer sites"


class Disposition(str, Enum):
    SEQUENCING = "sequencing"
    DISCARD = "discard"
    QC = "QC"
    STORED = "stored intermediate"
    UNRESOLVED = "unresolved"


class AnnotationSlot(str, Enum):
    """Where a graph-owned explanation appears on a protocol page."""
    INTRO = "intro"
    BEFORE_WORKFLOW = "before_workflow"
    AFTER_WORKFLOW = "after_workflow"
    SEQUENCING = "sequencing"


class DetailStyle(str, Enum):
    PANEL = "panel"
    FOLDED = "folded"


_ID = re.compile(r"^[a-z][a-z0-9_]*$")
_CONTRIBUTING = {
    InputRole.TEMPLATE, InputRole.STRUCTURAL, InputRole.ANALYTE, InputRole.SPIKE_IN,
}


def _identifier(value: str, what: str) -> str:
    if not _ID.fullmatch(value):
        raise ValueError(f"{what} must match {_ID.pattern!r}, got {value!r}")
    return value


@dataclass(frozen=True)
class Material:
    id: str
    label: str
    kind: MaterialKind
    state: Any = None
    primary: bool = False
    source: str = ""
    inferred: bool = False
    chemical_state: Any = None

    def __post_init__(self) -> None:
        _identifier(self.id, "material id")
        if not self.label.strip():
            raise ValueError("material needs a label")
        if self.primary and self.kind in {
                MaterialKind.CATALYST, MaterialKind.SELECTION, MaterialKind.CONDITION}:
            raise ValueError("a process reagent cannot be the primary molecular input")
        if self.inferred and not self.source.strip():
            raise ValueError("an inferred material needs a source-scoped explanation")


@dataclass(frozen=True)
class Input:
    material: str
    role: InputRole

    def __post_init__(self) -> None:
        _identifier(self.material, "input material id")


@dataclass(frozen=True)
class Step:
    id: str
    label: str
    kind: StepKind
    inputs: tuple[Input, ...]
    outputs: tuple[str, ...]
    geometry: Any = None
    source: str = ""
    changes_molecule: bool = True
    reason: str = ""
    relations: tuple[MolecularRelation, ...] = ()

    def __post_init__(self) -> None:
        _identifier(self.id, "step id")
        if not self.label.strip():
            raise ValueError("step needs a label")
        if self.kind is StepKind.MAGIC and not self.reason.strip():
            raise ValueError(f"magic step {self.id!r} needs a reason")
        if self.kind is StepKind.MAGIC and self.relations:
            raise ValueError(f"magic step {self.id!r} cannot claim checked relations")
        if len(set(self.relations)) != len(self.relations):
            raise ValueError(f"step {self.id!r} repeats a molecular relation")
        if any(not isinstance(relation, MolecularRelation) for relation in self.relations):
            raise ValueError(f"step {self.id!r} has an unknown molecular relation")
        if not self.inputs:
            raise ValueError(f"step {self.id!r} needs an input")
        if not self.outputs:
            raise ValueError(f"step {self.id!r} needs an output")
        if len(set(self.outputs)) != len(self.outputs):
            raise ValueError(f"step {self.id!r} repeats an output")
        for output in self.outputs:
            _identifier(output, "output material id")


@dataclass(frozen=True)
class Terminal:
    material: str
    disposition: Disposition
    reason: str = ""

    def __post_init__(self) -> None:
        _identifier(self.material, "terminal material id")
        if self.disposition in {Disposition.DISCARD, Disposition.STORED,
                                Disposition.UNRESOLVED} and not self.reason.strip():
            raise ValueError(f"a {self.disposition.value} terminal needs an explanation")


@dataclass(frozen=True)
class SequencingPlan:
    """How one terminal library is acquired.

    Exact plans carry the final ``Construct`` and sequencing primers.  If an authority
    withholds that information, ``limitation`` is mandatory and the endpoint remains
    explicit without inventing primer sites.
    """

    library: Any = None
    primers: tuple[Any, ...] = ()
    acquired_roles: tuple[str, ...] = ()
    read_lengths: Mapping[str, int] = field(default_factory=dict)
    read_uses: Mapping[str, str] = field(default_factory=dict)
    limitation: str = ""

    def __post_init__(self) -> None:
        exact = self.library is not None and bool(self.primers)
        if not exact and not self.limitation.strip():
            raise ValueError("sequencing plan needs exact primer geometry or a limitation")
        if self.primers and self.library is None:
            raise ValueError("sequencing primers need a final library")
        if len(set(self.acquired_roles)) != len(self.acquired_roles):
            raise ValueError("acquired sequencing roles must be unique")


@dataclass(frozen=True)
class SourceItem:
    id: str
    label: str
    source: str
    required: bool = True

    def __post_init__(self) -> None:
        _identifier(self.id, "source-item id")
        if not self.label.strip() or not self.source.strip():
            raise ValueError("source item needs a label and authoritative source")


@dataclass(frozen=True)
class SourceCoverage:
    item: str
    disposition: str
    target: str = ""
    reason: str = ""

    def __post_init__(self) -> None:
        _identifier(self.item, "covered source-item id")
        if self.disposition not in {"modelled", "not applicable", "excluded"}:
            raise ValueError("source coverage must be modelled, not applicable or excluded")
        if self.disposition == "modelled" and not self.target:
            raise ValueError("modelled source coverage needs a graph target")
        if self.disposition != "modelled" and not self.reason.strip():
            raise ValueError("excluded source coverage needs a reason")


@dataclass(frozen=True)
class Highlight:
    target: str
    title: str
    explanation: str
    slot: AnnotationSlot = AnnotationSlot.AFTER_WORKFLOW


@dataclass(frozen=True)
class DetailView:
    target: str
    title: str
    caption: str
    state: Any = None  # defaults to the targeted material's drawing state
    slot: AnnotationSlot = AnnotationSlot.AFTER_WORKFLOW
    style: DetailStyle = DetailStyle.PANEL


@dataclass(frozen=True)
class SourceClaim:
    location: str
    claim: str


@dataclass(frozen=True)
class SourceDiscrepancy:
    target: str
    title: str
    claims: tuple[SourceClaim, ...]
    model_choice: str
    evidence: str
    uncertainty: str
    slot: AnnotationSlot = AnnotationSlot.INTRO


GraphAnnotation = Highlight | DetailView | SourceDiscrepancy


class ProtocolGraph:
    """A typed protocol hypergraph built strictly in topological order."""

    def __init__(self, name: str, *, authority: str = ""):
        if not name.strip():
            raise ValueError("protocol graph needs a name")
        self.name = name
        self.authority = authority
        self.materials: dict[str, Material] = {}
        self.steps: dict[str, Step] = {}
        self.terminals: dict[str, Terminal] = {}
        self.sequencing: dict[str, SequencingPlan] = {}
        self.source_items: dict[str, SourceItem] = {}
        self.source_coverage: dict[str, SourceCoverage] = {}
        self.annotations: list[GraphAnnotation] = []
        self._producer: dict[str, str] = {}
        self._consumers: dict[str, list[str]] = {}

    def add_input(self, material: Material) -> Material:
        """Declare a source material or reagent with no producer inside the protocol."""
        self._add_material(material)
        return material

    def input(self, identifier: str, label: str, kind: MaterialKind, *,
              state: Any = None, chemical_state: Any = None,
              primary: bool = False, source: str | None = None,
              inferred: bool = False) -> Material:
        """Construct and declare an external material using the graph authority by default."""
        material = Material(identifier, label, kind, state=state,
                            chemical_state=chemical_state, primary=primary,
                            source=self.authority if source is None else source,
                            inferred=inferred)
        return self.add_input(material)

    def product(self, identifier: str, label: str, kind: MaterialKind, *,
                state: Any = None, chemical_state: Any = None,
                source: str | None = None,
                inferred: bool = False) -> Material:
        """Construct a prospective step output; :meth:`react` installs its producer."""
        return Material(identifier, label, kind, state=state,
                        chemical_state=chemical_state,
                        source=self.authority if source is None else source,
                        inferred=inferred)

    def react(self, identifier: str, label: str, kind: StepKind,
              inputs: Iterable[Input | tuple[Material | str, InputRole]],
              outputs: Iterable[Material], *, geometry: Any = None,
              source: str | None = None, changes_molecule: bool = True,
              reason: str = "", relations: Iterable[MolecularRelation] = ()) -> Step:
        """Construct a typed step while preserving explicit role-labelled ports."""
        ports = tuple(
            value if isinstance(value, Input)
            else Input(value[0].id if isinstance(value[0], Material) else value[0], value[1])
            for value in inputs)
        made = tuple(outputs)
        step = Step(identifier, label, kind, ports, tuple(x.id for x in made),
                    geometry=geometry,
                    source=self.authority if source is None else source,
                    changes_molecule=changes_molecule, reason=reason,
                    relations=tuple(relations))
        return self.add_step(step, made)

    def sequence(self, material: Material | str, plan: SequencingPlan) -> Terminal:
        """Declare a sequencing terminal without repeating its material identifier."""
        identifier = material.id if isinstance(material, Material) else material
        return self.mark_terminal(Terminal(identifier, Disposition.SEQUENCING), plan)

    def discard(self, material: Material | str, reason: str) -> Terminal:
        """Declare a discarded terminal; the validator still requires waste typing."""
        identifier = material.id if isinstance(material, Material) else material
        return self.mark_terminal(Terminal(identifier, Disposition.DISCARD, reason))

    def _add_material(self, material: Material) -> None:
        if material.id in self.materials:
            raise ValueError(f"duplicate material id {material.id!r}")
        if material.id in self.steps:
            raise ValueError(f"material id {material.id!r} collides with a step")
        self.materials[material.id] = material
        self._consumers[material.id] = []

    def add_step(self, step: Step, outputs: Iterable[Material]) -> Step:
        """Add a reaction and its newly created outputs.

        Inputs must already exist and outputs must be new.  This topological construction
        rule means a cycle or a material with two producers cannot be represented.
        """
        if step.id in self.steps or step.id in self.materials:
            raise ValueError(f"duplicate or colliding step id {step.id!r}")
        missing = [x.material for x in step.inputs if x.material not in self.materials]
        if missing:
            raise ValueError(f"step {step.id!r} refers to unknown inputs: {', '.join(missing)}")
        allowed = {
            InputRole.TEMPLATE: {MaterialKind.BIOLOGICAL, MaterialKind.MOLECULAR,
                                 MaterialKind.OLIGO, MaterialKind.ADAPTER,
                                 MaterialKind.LIBRARY, MaterialKind.QC},
            InputRole.STRUCTURAL: {MaterialKind.OLIGO, MaterialKind.ADAPTER,
                                   MaterialKind.MOLECULAR},
            InputRole.CATALYST: {MaterialKind.CATALYST},
            InputRole.SELECTION: {MaterialKind.SELECTION},
            InputRole.ANALYTE: {MaterialKind.BIOLOGICAL, MaterialKind.MOLECULAR},
            InputRole.CONDITION: {MaterialKind.CONDITION},
            InputRole.SPIKE_IN: {MaterialKind.BIOLOGICAL, MaterialKind.MOLECULAR,
                                 MaterialKind.OLIGO},
            InputRole.SELF_PRIMER: {MaterialKind.MOLECULAR,MaterialKind.OLIGO},
        }
        for port in step.inputs:
            kind = self.materials[port.material].kind
            if kind not in allowed[port.role]:
                raise ValueError(
                    f"{port.material!r} is {kind.value}, invalid for {port.role.value}")
        made = tuple(outputs)
        if tuple(x.id for x in made) != step.outputs:
            raise ValueError("step output ids and supplied material objects differ")
        for material in made:
            if material.id in self.materials or material.id in self.steps:
                raise ValueError(f"output {material.id!r} already exists")
        self._validate_step_shape(step)
        self.steps[step.id] = step
        for port in step.inputs:
            self._consumers[port.material].append(step.id)
        for material in made:
            self._add_material(material)
            self._producer[material.id] = step.id
        return step

    def mark_terminal(self, terminal: Terminal,
                      sequencing: SequencingPlan | None = None) -> Terminal:
        if terminal.material not in self.materials:
            raise ValueError(f"terminal refers to unknown material {terminal.material!r}")
        if terminal.material in self.terminals:
            raise ValueError(f"material {terminal.material!r} already has a disposition")
        if (terminal.disposition is Disposition.SEQUENCING) != (sequencing is not None):
            raise ValueError("a sequencing terminal needs exactly one sequencing plan")
        self.terminals[terminal.material] = terminal
        if sequencing is not None:
            self.sequencing[terminal.material] = sequencing
        return terminal

    def add_source_item(self, item: SourceItem) -> SourceItem:
        if item.id in self.source_items:
            raise ValueError(f"duplicate source item {item.id!r}")
        self.source_items[item.id] = item
        return item

    def annotate(self, annotation: GraphAnnotation) -> GraphAnnotation:
        """Attach page emphasis to an existing node without changing chemistry claims."""
        if not isinstance(annotation, (Highlight, DetailView, SourceDiscrepancy)):
            raise TypeError("unknown graph annotation")
        if annotation.target not in self.materials and annotation.target not in self.steps:
            raise ValueError(f"annotation names unknown graph target {annotation.target!r}")
        if not annotation.title.strip():
            raise ValueError("annotation needs a title")
        if not isinstance(annotation.slot, AnnotationSlot):
            raise ValueError("annotation needs a known page slot")
        if isinstance(annotation, Highlight) and not annotation.explanation.strip():
            raise ValueError("highlight needs an explanation")
        if isinstance(annotation, DetailView):
            if not annotation.caption.strip() or not isinstance(annotation.style, DetailStyle):
                raise ValueError("detail view needs a caption and known style")
            if annotation.state is None and (annotation.target not in self.materials or
                                             self.materials[annotation.target].state is None):
                raise ValueError("detail view needs a drawing state")
            if annotation.style is DetailStyle.FOLDED:
                state = annotation.state or self.materials[annotation.target].state
                if not callable(getattr(state, "folded_panel", None)):
                    raise ValueError("folded detail view needs folded_panel(caption)")
        if isinstance(annotation, SourceDiscrepancy):
            if len(annotation.claims) < 2 or any(
                    not claim.location.strip() or not claim.claim.strip()
                    for claim in annotation.claims):
                raise ValueError("source discrepancy needs at least two located claims")
            if not all((annotation.model_choice.strip(), annotation.evidence.strip(),
                        annotation.uncertainty.strip())):
                raise ValueError("source discrepancy needs choice, evidence and uncertainty")
        self.annotations.append(annotation)
        return annotation

    def cover_source_item(self, coverage: SourceCoverage) -> SourceCoverage:
        if coverage.item not in self.source_items:
            raise ValueError(f"coverage names undeclared source item {coverage.item!r}")
        if coverage.item in self.source_coverage:
            raise ValueError(f"source item {coverage.item!r} is covered twice")
        if coverage.disposition == "modelled":
            targets = set(self.materials) | set(self.steps)
            if coverage.target not in targets:
                raise ValueError(f"coverage target {coverage.target!r} is not in the graph")
        self.source_coverage[coverage.item] = coverage
        return coverage

    def declare_reviewed_scope(self, targets: Iterable[str] | Mapping[str, str],
                               source: str) -> None:
        """Declare the independently reviewed source inventory for graph targets.

        ``targets`` remains an explicit curator-written list: deriving it from the graph
        would only prove that the graph agrees with itself.  This helper centralizes the
        repetitive one-source-item/one-target registration without changing that review
        boundary.
        """
        if not source.strip():
            raise ValueError("reviewed scope needs an authoritative source")
        pairs = (targets.items() if isinstance(targets, Mapping)
                 else (("source_" + target, target) for target in targets))
        seen_items: set[str] = set()
        valid = set(self.materials) | set(self.steps)
        for item_id, target in pairs:
            if item_id in seen_items:
                raise ValueError(f"reviewed scope repeats source item {item_id!r}")
            if target not in valid:
                raise ValueError(f"reviewed scope target {target!r} is not in the graph")
            seen_items.add(item_id)
            self.add_source_item(SourceItem(
                item_id, target.replace("_", " "), source))
            self.cover_source_item(SourceCoverage(item_id, "modelled", target))

    def _validate_step_shape(self, step: Step) -> None:
        roles = [p.role for p in step.inputs]
        n = len(step.outputs)

        def need(role: InputRole, message: str) -> None:
            if role not in roles:
                raise ValueError(f"{step.kind.value} step {step.id!r} {message}")

        if step.kind in {StepKind.PCR, StepKind.LIGATION,
                         StepKind.CIRCULARIZATION} and step.geometry is None:
            raise ValueError(f"{step.kind.value} step {step.id!r} needs structured geometry")
        if step.kind is StepKind.MAGIC:
            need(InputRole.TEMPLATE, "needs a substrate")
        elif step.kind is StepKind.PCR:
            need(InputRole.TEMPLATE, "needs a template")
            primer_inputs = roles.count(InputRole.STRUCTURAL)
            symmetric = getattr(step.geometry,"symmetric_species",False)
            if primer_inputs < 1 or (primer_inputs < 2 and not symmetric):
                raise ValueError(
                    f"PCR step {step.id!r} needs both primer inputs, or one exact "
                    "primer species with symmetric inward sites")
            if getattr(step.geometry, "availability", None) not in {
                    "exact", "regional", "unavailable"}:
                raise ValueError(f"PCR step {step.id!r} has invalid geometry")
        elif step.kind is StepKind.LIGATION:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.STRUCTURAL, "needs a ligating reagent")
            if getattr(step.geometry, "availability", None) not in {"known", "unavailable"}:
                raise ValueError(f"ligation step {step.id!r} has invalid geometry")
        elif step.kind is StepKind.PROXIMITY_LIGATION:
            need(InputRole.TEMPLATE, "needs restriction-fragment substrates")
            need(InputRole.CATALYST, "needs a ligase")
        elif step.kind is StepKind.DIGESTION:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.CATALYST, "needs a nuclease")
        elif step.kind is StepKind.TAGMENTATION:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.STRUCTURAL, "needs a loaded adapter")
            need(InputRole.CATALYST, "needs a transposase")
            if n < 2:
                raise ValueError(
                    "tagmentation must expose library-eligible and non-library outputs")
        elif step.kind is StepKind.ADAPTER_ATTACHMENT:
            need(InputRole.TEMPLATE, "needs a tagged molecular substrate")
            need(InputRole.STRUCTURAL, "needs an attaching adapter")
        elif step.kind is StepKind.REVERSE_TRANSCRIPTION:
            need(InputRole.TEMPLATE, "needs an RNA template")
            if not ({InputRole.STRUCTURAL,InputRole.SELF_PRIMER} & set(roles)):
                raise ValueError(
                    f"reverse transcription step {step.id!r} needs an external or intramolecular RT primer")
            need(InputRole.CATALYST, "needs reverse transcriptase")
        elif step.kind is StepKind.TEMPLATE_SWITCHING:
            need(InputRole.TEMPLATE, "needs a nascent cDNA substrate")
            need(InputRole.STRUCTURAL, "needs a template-switch oligo")
            need(InputRole.CATALYST, "needs template-switching reverse transcriptase")
        elif step.kind is StepKind.TRANSCRIPTION:
            need(InputRole.TEMPLATE, "needs a DNA template")
            need(InputRole.CATALYST, "needs an RNA polymerase")
        elif step.kind is StepKind.RUN_ON:
            need(InputRole.TEMPLATE, "needs engaged transcription complexes")
            need(InputRole.STRUCTURAL, "needs labelled run-on nucleotide")
            need(InputRole.CATALYST, "needs engaged RNA polymerase")
        elif step.kind is StepKind.HYBRIDIZATION:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.STRUCTURAL, "needs a hybridizing reagent")
        elif step.kind is StepKind.BINDING:
            need(InputRole.TEMPLATE, "needs a binding substrate")
            if not ({InputRole.SELECTION, InputRole.ANALYTE} & set(roles)):
                raise ValueError(
                    f"binding step {step.id!r} needs a binding probe or measured analyte")
            if n < 2:
                raise ValueError("binding must expose bound and unbound outcomes")
        elif step.kind is StepKind.FOLDING:
            need(InputRole.TEMPLATE,"needs a molecule to fold")
            need(InputRole.CONDITION,"needs a folding condition")
        elif step.kind in {StepKind.EXTENSION, StepKind.TAILING}:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.CATALYST, "needs an enzyme")
        elif step.kind in {StepKind.END_REPAIR, StepKind.PHOSPHORYLATION,
                           StepKind.DEPHOSPHORYLATION}:
            need(InputRole.TEMPLATE, "needs a substrate")
            need(InputRole.CATALYST, "needs an enzyme")
        elif step.kind is StepKind.NICK_SEALING:
            need(InputRole.TEMPLATE, "needs a nicked substrate")
            need(InputRole.CATALYST, "needs a ligase")
        elif step.kind is StepKind.CHEMICAL_CONVERSION:
            need(InputRole.TEMPLATE, "needs a substrate")
            if not ({InputRole.CATALYST, InputRole.CONDITION} & set(roles)):
                raise ValueError(
                    f"chemical conversion step {step.id!r} needs a catalyst or condition")
        elif step.kind is StepKind.CROSS_LINKING:
            need(InputRole.TEMPLATE, "needs material to cross-link")
            need(InputRole.CONDITION, "needs a cross-linking reagent")
        elif step.kind is StepKind.PERMEABILIZATION:
            need(InputRole.TEMPLATE, "needs cells or nuclei")
            need(InputRole.CONDITION, "needs a permeabilization condition")
        elif step.kind is StepKind.SOLUBILIZATION:
            need(InputRole.TEMPLATE, "needs material to solubilize")
            need(InputRole.CONDITION, "needs a solubilization condition")
            if n < 2:
                raise ValueError("solubilization must expose soluble and insoluble outcomes")
        elif step.kind is StepKind.QUENCHING:
            need(InputRole.TEMPLATE, "needs an active reaction")
            need(InputRole.CONDITION, "needs a quenching condition")
        elif step.kind is StepKind.THERMAL_TREATMENT:
            need(InputRole.TEMPLATE, "needs material to heat or cool")
            need(InputRole.CONDITION, "needs a thermal condition")
        elif step.kind is StepKind.FRAGMENTATION:
            need(InputRole.TEMPLATE, "needs material to fragment")
            if not ({InputRole.CATALYST, InputRole.CONDITION} & set(roles)):
                raise ValueError(
                    f"fragmentation step {step.id!r} needs an enzyme or physical condition")
        elif step.kind is StepKind.DEPROTEINIZATION:
            need(InputRole.TEMPLATE, "needs a protein-bearing substrate")
            if not ({InputRole.CATALYST, InputRole.CONDITION} & set(roles)):
                raise ValueError(
                    f"deproteinization step {step.id!r} needs an enzyme or condition")
        elif step.kind is StepKind.EXTRACTION:
            need(InputRole.TEMPLATE, "needs material to extract")
            if not ({InputRole.CONDITION, InputRole.SELECTION} & set(roles)):
                raise ValueError(
                    f"extraction step {step.id!r} needs an extraction condition or selection reagent")
        elif step.kind is StepKind.STRAND_SEPARATION:
            need(InputRole.TEMPLATE, "needs a duplex substrate")
            need(InputRole.CONDITION, "needs a denaturing condition")
            if n < 2:
                raise ValueError("strand separation must expose both strand outcomes")
        elif step.kind is StepKind.INTEGRATION:
            need(InputRole.TEMPLATE, "needs an integrating construct")
            need(InputRole.CATALYST, "needs an integration catalyst")
        elif step.kind is StepKind.GENOME_EDITING:
            need(InputRole.TEMPLATE, "needs a genomic target")
            need(InputRole.STRUCTURAL, "needs a guide reagent")
            need(InputRole.CATALYST, "needs a genome-editing catalyst")
        elif step.kind is StepKind.CIRCULARIZATION:
            need(InputRole.TEMPLATE, "needs a linear substrate")
            need(InputRole.CATALYST, "needs a ligase")
            if getattr(step.geometry, "availability", None) not in {"known", "unavailable"}:
                raise ValueError(f"circularization step {step.id!r} has invalid geometry")
        elif step.kind in {StepKind.SELECTION, StepKind.AFFINITY_CAPTURE}:
            need(InputRole.TEMPLATE, "needs material to select")
            need(InputRole.SELECTION, "needs a selection reagent")
            if n < 2:
                raise ValueError("selection must expose retained and rejected outputs")
        elif step.kind is StepKind.SAMPLING:
            need(InputRole.TEMPLATE, "needs material to sample")
            if not ({InputRole.CONDITION, InputRole.SELECTION} & set(roles)):
                raise ValueError("sampling needs a bottleneck condition or selection rule")
            if n < 2:
                raise ValueError("sampling must expose retained and excluded material")
        elif step.kind in {StepKind.SPLIT, StepKind.PARTITION}:
            need(InputRole.TEMPLATE, "needs material to divide")
            if n < 2:
                raise ValueError(f"{step.kind.value} needs at least two outputs")
        elif step.kind is StepKind.POOL:
            if roles.count(InputRole.TEMPLATE) < 2 or n != 1:
                raise ValueError("pool needs at least two templates and exactly one output")
        elif step.kind is StepKind.CLEANUP:
            need(InputRole.TEMPLATE, "needs material to clean")
            if n < 2:
                raise ValueError("cleanup must expose retained and discarded outputs")

    def _downstream_terminals(self, material_id: str) -> set[str]:
        found: set[str] = set()
        stack = [material_id]
        seen: set[str] = set()
        while stack:
            material = stack.pop()
            if material in seen:
                continue
            seen.add(material)
            if material in self.terminals:
                found.add(material)
            for step_id in self._consumers[material]:
                stack.extend(self.steps[step_id].outputs)
        return found

    def producer_of(self, material_id: str) -> Step | None:
        """The unique producing step, or ``None`` for a declared protocol input."""
        if material_id not in self.materials:
            raise KeyError(material_id)
        step_id = self._producer.get(material_id)
        return self.steps[step_id] if step_id is not None else None

    def consumers_of(self, material_id: str) -> tuple[Step, ...]:
        """Every step that consumes a material, in construction order."""
        if material_id not in self.materials:
            raise KeyError(material_id)
        return tuple(self.steps[x] for x in self._consumers[material_id])

    def _contributing_origins(self, material_id: str) -> set[str]:
        origins: set[str] = set()
        stack = [material_id]
        seen: set[str] = set()
        while stack:
            material = stack.pop()
            if material in seen:
                continue
            seen.add(material)
            producer = self._producer.get(material)
            if producer is None:
                origins.add(material)
                continue
            step = self.steps[producer]
            stack.extend(p.material for p in step.inputs if p.role in _CONTRIBUTING)
        return origins

    def publication_problems(self) -> list[str]:
        """Return every reason this graph is not eligible for a finished public page."""
        problems: list[str] = []
        if not self.authority.strip():
            problems.append("no defining paper or commercial authority")
        sequencing = [t.material for t in self.terminals.values()
                      if t.disposition is Disposition.SEQUENCING]
        if not sequencing:
            problems.append("no terminal output is designated for sequencing")
        if not self.source_items:
            problems.append("no authoritative source items declare the reviewed scope")

        for material_id, material in self.materials.items():
            consumed = bool(self._consumers[material_id])
            terminal = material_id in self.terminals
            if not consumed and not terminal:
                problems.append(f"material {material_id!r} is neither consumed nor terminal")
            if consumed and terminal:
                problems.append(f"terminal material {material_id!r} is consumed again")
            if material.primary and material_id in self._producer:
                problems.append(f"primary input {material_id!r} has an internal producer")
            if not material.source.strip():
                problems.append(f"material {material_id!r} has no source attribution")

        for terminal in self.terminals.values():
            kind = self.materials[terminal.material].kind
            if terminal.disposition is Disposition.DISCARD and kind is not MaterialKind.WASTE:
                problems.append(f"discard terminal {terminal.material!r} is not typed as waste")
            if terminal.disposition is Disposition.QC and kind is not MaterialKind.QC:
                problems.append(f"QC terminal {terminal.material!r} is not typed as QC material")

        for step in self.steps.values():
            if not any(self._downstream_terminals(out) for out in step.outputs):
                problems.append(f"step {step.id!r} does not reach a terminal outcome")
            if step.kind is StepKind.GENERIC and step.changes_molecule:
                problems.append(f"construct-changing step {step.id!r} has generic chemistry")
            if not step.source.strip():
                problems.append(f"step {step.id!r} has no source attribution")

        for material_id in sequencing:
            material = self.materials[material_id]
            if material.kind is not MaterialKind.LIBRARY:
                problems.append(f"sequencing output {material_id!r} is not a library")
            origins = self._contributing_origins(material_id)
            if not any(self.materials[x].primary for x in origins):
                problems.append(f"sequencing output {material_id!r} lacks primary-template provenance")
            plan = self.sequencing[material_id]
            if plan.library is not None and plan.primers:
                # Local import avoids a module cycle: seqprimers imports chemdraw, whose
                # Workflow in turn uses ProtocolGraph.
                import seqprimers
                required = plan.acquired_roles or tuple(p.role for p in plan.primers)
                for error in seqprimers.verify(plan.library, plan.primers, required):
                    problems.append(f"sequencing output {material_id!r}: {error}")
                try:
                    seqprimers.read_contents(plan.library, plan.primers,
                                             dict(plan.read_lengths), dict(plan.read_uses),
                                             acquired_roles=required)
                except ValueError as exc:
                    problems.append(f"sequencing output {material_id!r}: {exc}")

        for material in self.materials.values():
            if material.kind in {MaterialKind.OLIGO, MaterialKind.ADAPTER}:
                terminals = self._downstream_terminals(material.id)
                if not terminals:
                    problems.append(f"structural input {material.id!r} has no terminal fate")

        for item in self.source_items.values():
            if item.required and item.id not in self.source_coverage:
                problems.append(f"authoritative source item {item.id!r} is not accounted for")
        covered_targets = {coverage.target for coverage in self.source_coverage.values()
                           if coverage.disposition == "modelled"}
        required_targets = (set(self.steps)
                            | {key for key in self.materials if key not in self._producer}
                            | set(sequencing))
        for target in sorted(required_targets - covered_targets):
            problems.append(f"graph input/step {target!r} has no source-inventory item")

        # Weak connectivity is checked across material co-participation in each hyperedge.
        # This prevents an independently valid but unrelated discard/QC mini-protocol from
        # being smuggled into the same named protocol.
        if self.materials:
            neighbours = {key: set() for key in self.materials}
            for step in self.steps.values():
                participants = [x.material for x in step.inputs] + list(step.outputs)
                for material in participants:
                    neighbours[material].update(x for x in participants if x != material)
            seen: set[str] = set()
            stack = [next(iter(self.materials))]
            while stack:
                material = stack.pop()
                if material in seen:
                    continue
                seen.add(material)
                stack.extend(neighbours[material] - seen)
            if len(seen) != len(self.materials):
                missing = ", ".join(sorted(set(self.materials) - seen))
                problems.append(f"protocol graph is disconnected: {missing}")
        return problems

    def require_publishable(self) -> "ProtocolGraph":
        problems = self.publication_problems()
        if problems:
            raise ValueError(f"{self.name} is incomplete:\n- " + "\n- ".join(problems))
        return self

    @property
    def publishable(self) -> bool:
        return not self.publication_problems()
