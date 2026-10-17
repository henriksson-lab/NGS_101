"""Bounded chemistry and reaction operations for RNA end/readout protocols."""
from __future__ import annotations

from dataclasses import dataclass

from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind
from protocol_model import molecule, role


@dataclass(frozen=True)
class StopReadout:
    reagent: str
    target: str
    control: str

    def __post_init__(self):
        if not self.reagent or not self.target or not self.control:
            raise ValueError("a modification-stop readout requires reagent, target and control")


@dataclass(frozen=True)
class BaseConversionReadout:
    reagent: str
    converted_base: str
    product: str
    read_as: str
    protected_base: str

    def __post_init__(self):
        if self.converted_base == self.protected_base:
            raise ValueError("converted and protected bases must be distinct chemical states")
        if not all((self.reagent, self.product, self.read_as)):
            raise ValueError("base-conversion chemistry must declare reagent, product and readout")


@dataclass(frozen=True)
class SelectiveCleavageReadout:
    target_modifications: tuple[str, ...]
    cleavage_reagent: str
    selected_end: str
    signal_offset: int

    def __post_init__(self):
        if self.selected_end != "5-prime phosphate":
            raise ValueError("adapter-positive cleavage readout must select a newly made 5-prime phosphate")
        if self.signal_offset != 1 or not self.target_modifications:
            raise ValueError("cleavage signal must be the N+1 fragment start for declared targets")


def _finish_pcr(graph, prefix, template, final, endpoint, reviewed, *, label,
                geometry=None, primers=(), acquired_roles=(), read_lengths=None,
                read_uses=None, step_label=None, forward_label=None,
                reverse_label=None, polymerase_label=None, forward_state=None,
                reverse_state=None, primer_source=None, primers_inferred=False,
                unavailable_reason=None):
    """Amplify ``template`` into the caller's ``final`` construct and sequence it.

    The default wording names roles only: it does not claim that the PCR adds a
    sample index or which primers and polymerase were used.  Callers name source
    reagents, pass printed primer states (``forward_state``/``reverse_state``) and
    either an exact ``geometry`` or the source-specific ``unavailable_reason``.
    ``primers_inferred`` marks primer inputs whose identity is not source-printed.
    """
    forward = graph.input(prefix + "_forward", forward_label or label + " forward PCR primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source,
                          inferred=primers_inferred)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " reverse PCR primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source,
                          inferred=primers_inferred)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", step_label or "Amplify the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=(geometry or PCRGeometry.unavailable(
                    final, reason=unavailable_reason or "The source establishes the readout roles but not one universal source-exact amplicon geometry.")))
    if primers:
        graph.sequence(library, SequencingPlan(
            library=final, primers=tuple(primers),
            acquired_roles=tuple(acquired_roles),
            read_lengths=read_lengths or {}, read_uses=read_uses or {}))
    else:
        graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


_REAGENT_ROLE = {
    StepKind.LIGATION: ("_adapter", " adapter or ligating oligo", MaterialKind.ADAPTER,
                        InputRole.STRUCTURAL),
    StepKind.REVERSE_TRANSCRIPTION: ("_primer", " RT primer", MaterialKind.OLIGO,
                                     InputRole.STRUCTURAL),
    StepKind.TEMPLATE_SWITCHING: ("_oligo", " template oligo", MaterialKind.OLIGO,
                                  InputRole.STRUCTURAL),
    StepKind.FRAGMENTATION: ("_condition", " fragmentation condition", MaterialKind.CONDITION,
                             InputRole.CONDITION),
    StepKind.CHEMICAL_CONVERSION: ("_condition", " condition", MaterialKind.CONDITION,
                                   InputRole.CONDITION),
    StepKind.THERMAL_TREATMENT: ("_condition", " condition", MaterialKind.CONDITION,
                                 InputRole.CONDITION),
    StepKind.DEPROTEINIZATION: ("_condition", " condition", MaterialKind.CONDITION,
                                InputRole.CONDITION),
    StepKind.SELECTION: ("_selection", " selection criterion", MaterialKind.SELECTION,
                         InputRole.SELECTION),
}

_CATALYST_ROLE = {
    StepKind.LIGATION: ("_ligase", " ligase"),
    StepKind.REVERSE_TRANSCRIPTION: ("_enzyme", " reverse transcriptase"),
    StepKind.TEMPLATE_SWITCHING: ("_enzyme", " template-switching polymerase"),
    **{kind: ("_enzyme", " enzyme") for kind in (
        StepKind.DIGESTION, StepKind.END_REPAIR, StepKind.PHOSPHORYLATION,
        StepKind.DEPHOSPHORYLATION, StepKind.CIRCULARIZATION, StepKind.TAILING)},
}


def _apply_linear_step(graph, current, identifier, label, kind, reviewed, *,
                       reagent_label=None, reagent_state=None, reagent_source=None,
                       catalyst_label=None, product_state=None, geometry=None,
                       discarded_label=None):
    """Apply one role-level step to ``current`` and return the retained product.

    Optional keywords let a caller supply source facts without changing the step
    identifiers: the non-catalyst reagent's label, exact state and source, the
    catalyst label, the product state and exact geometry.  SELECTION and CLEANUP
    steps keep the input molecule unchanged and discard the excluded material.
    """
    inputs = [(current, InputRole.TEMPLATE)]
    if kind not in _REAGENT_ROLE and kind not in _CATALYST_ROLE and kind is not StepKind.CLEANUP:
        raise ValueError(f"unsupported linear RNA-end step {kind}")
    if kind in _REAGENT_ROLE:
        suffix, text, material, port = _REAGENT_ROLE[kind]
        reagent = graph.input(identifier + suffix, reagent_label or label + text, material,
                              state=reagent_state, source=reagent_source)
        inputs.append((reagent, port))
        reviewed.append(identifier + suffix)
    elif reagent_label or reagent_state is not None:
        raise ValueError(f"{kind.value} step {identifier!r} takes no structural reagent")
    if kind in _CATALYST_ROLE:
        suffix, text = _CATALYST_ROLE[kind]
        enzyme = graph.input(identifier + suffix, catalyst_label or label + text,
                             MaterialKind.CATALYST)
        inputs.append((enzyme, InputRole.CATALYST))
        reviewed.append(identifier + suffix)
    elif catalyst_label:
        raise ValueError(f"{kind.value} step {identifier!r} takes no catalyst")
    if geometry is None and kind is StepKind.LIGATION:
        geometry = AdapterLigationGeometry.unavailable(
            "The source defines the ligated molecular role; the role-level shared page does not assert one universal end geometry.")
    elif geometry is None and kind is StepKind.CIRCULARIZATION:
        geometry = CircularizationGeometry.unavailable(
            "The source establishes cDNA circularization; this shared role-level model does not assert one universal closing-end geometry.")
    if kind in {StepKind.SELECTION, StepKind.CLEANUP}:
        product = graph.product(identifier + "_product", label + " retained material",
                                MaterialKind.MOLECULAR,
                                state=current.state if product_state is None else product_state)
        excluded = graph.product(identifier + "_excluded",
                                 discarded_label or label + " excluded material",
                                 MaterialKind.WASTE)
        graph.react(identifier, label, kind, inputs, (product, excluded),
                    geometry=geometry, changes_molecule=False)
        graph.discard(excluded, "not carried into the retained library material")
        reviewed.append(identifier)
        return product
    product = graph.product(identifier + "_product", label + " product",
                            MaterialKind.MOLECULAR,
                            state=(product_state if product_state is not None else
                                   molecule(label + " product", (role("retained molecular product"),))))
    graph.react(identifier, label, kind, inputs, (product,), geometry=geometry)
    reviewed.append(identifier)
    return product
