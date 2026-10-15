"""Bounded scene, workflow, and library-PCR helpers for selection assays."""
from __future__ import annotations

from chemdraw import MolecularState, Scene, Segment, Workflow
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind


def _scene(label: str, *parts: Segment) -> Scene:
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    sc.labels("molecule")
    return sc

def _workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w

def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The source establishes the endpoint but withholds or varies the complete amplification-primer pair."))
    plan = (SequencingPlan(library=final, primers=tuple(primers)) if primers
            else SequencingPlan(limitation=endpoint))
    graph.sequence(library, plan)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library
