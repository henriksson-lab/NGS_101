"""Bounded PCR finishing operation for editing-outcome graphs."""
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                step_label=None, forward_label=None, reverse_label=None,
                polymerase_label=None, forward_state=None, reverse_state=None,
                primer_source=None, primers_inferred=False, geometry=None,
                unavailable_reason=None):
    """Endpoint PCR with role-level default labels.

    Callers name source primers and polymerase, pass printed primer states and an
    exact ``geometry`` or a source-specific ``unavailable_reason``.
    """
    forward = graph.input(prefix + "_forward", forward_label or label + " forward primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source,
                          inferred=primers_inferred)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " reverse primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source,
                          inferred=primers_inferred)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                            MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", step_label or label + " PCR", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry or PCRGeometry.unavailable(
                    final, reason=unavailable_reason or "Primer identities are represented in the final construct; the intermediate template is not an exact Construct."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                           limitation=""))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr",
                     prefix + "_library"))
    return library
