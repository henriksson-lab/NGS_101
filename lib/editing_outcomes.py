"""Bounded PCR finishing operation for editing-outcome graphs."""
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " forward primer", MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                            MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", label + " PCR", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="Primer identities are represented in the final construct; the intermediate template is not an exact Construct."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                           limitation=""))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr",
                     prefix + "_library"))
    return library
