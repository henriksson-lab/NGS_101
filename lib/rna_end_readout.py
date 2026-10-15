"""Bounded chemistry and reaction operations for RNA end/readout protocols."""
from __future__ import annotations

from dataclasses import dataclass

from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind
from protocol_model import molecule, role


class StopReadout:
    reagent: str
    target: str
    control: str

    def __post_init__(self):
        if not self.reagent or not self.target or not self.control:
            raise ValueError("a modification-stop readout requires reagent, target and control")


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
                read_uses=None):
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " indexed reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=(geometry or PCRGeometry.unavailable(
                    final, reason="The source establishes the readout roles but not one universal source-exact amplicon geometry.")))
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


def _apply_linear_step(graph, current, identifier, label, kind, reviewed):
    inputs = [(current, InputRole.TEMPLATE)]
    geometry = None
    if kind is StepKind.LIGATION:
        reagent = graph.input(identifier + "_adapter", label + " adapter or ligating oligo",
                              MaterialKind.ADAPTER)
        enzyme = graph.input(identifier + "_ligase", label + " ligase",
                             MaterialKind.CATALYST)
        inputs.extend(((reagent, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        geometry = AdapterLigationGeometry.unavailable(
            "The source defines the ligated molecular role; the role-level shared page does not assert one universal end geometry.")
        reviewed.extend((identifier + "_adapter", identifier + "_ligase"))
    elif kind is StepKind.REVERSE_TRANSCRIPTION:
        primer = graph.input(identifier + "_primer", label + " RT primer", MaterialKind.OLIGO)
        enzyme = graph.input(identifier + "_enzyme", label + " reverse transcriptase",
                             MaterialKind.CATALYST)
        inputs.extend(((primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        reviewed.extend((identifier + "_primer", identifier + "_enzyme"))
    elif kind is StepKind.TEMPLATE_SWITCHING:
        oligo = graph.input(identifier + "_oligo", label + " template oligo",
                            MaterialKind.OLIGO)
        enzyme = graph.input(identifier + "_enzyme", label + " template-switching polymerase",
                             MaterialKind.CATALYST)
        inputs.extend(((oligo, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)))
        reviewed.extend((identifier + "_oligo", identifier + "_enzyme"))
    elif kind in {StepKind.DIGESTION, StepKind.END_REPAIR, StepKind.PHOSPHORYLATION,
                  StepKind.DEPHOSPHORYLATION, StepKind.CIRCULARIZATION, StepKind.TAILING}:
        enzyme = graph.input(identifier + "_enzyme", label + " enzyme", MaterialKind.CATALYST)
        inputs.append((enzyme, InputRole.CATALYST))
        reviewed.append(identifier + "_enzyme")
        if kind is StepKind.CIRCULARIZATION:
            geometry = CircularizationGeometry.unavailable(
                "The source establishes cDNA circularization; this shared role-level model does not assert one universal closing-end geometry.")
    elif kind is StepKind.FRAGMENTATION:
        condition = graph.input(identifier + "_condition", label + " fragmentation condition",
                                MaterialKind.CONDITION)
        inputs.append((condition, InputRole.CONDITION))
        reviewed.append(identifier + "_condition")
    elif kind in {StepKind.CHEMICAL_CONVERSION, StepKind.THERMAL_TREATMENT,
                  StepKind.DEPROTEINIZATION}:
        condition = graph.input(identifier + "_condition", label + " condition",
                                MaterialKind.CONDITION)
        inputs.append((condition, InputRole.CONDITION))
        reviewed.append(identifier + "_condition")
    else:
        raise ValueError(f"unsupported linear RNA-end step {kind}")
    product = graph.product(identifier + "_product", label + " product",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " product", (role("retained molecular product"),)))
    graph.react(identifier, label, kind, inputs, (product,), geometry=geometry)
    reviewed.append(identifier)
    return product
