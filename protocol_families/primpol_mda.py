"""Shared TthPrimPol/phi29 reactions for the two TruePrime input workflows."""

from __future__ import annotations

from protocol_model import molecule, role
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind


PRIMER_LENGTH_NT = (7, 9)
AMPLIFICATION_C = 30
INACTIVATION_C = 65


def add_amplification(graph, template, *, authority: str, duration: str):
    """Prime exposed DNA de novo, then extend it by strand displacement."""
    primpol = graph.input("tth_primpol", "TthPrimPol", MaterialKind.CATALYST)
    primed = graph.product(
        "primed", "templates bearing de novo 7–9-nt DNA primers",
        MaterialKind.MOLECULAR,
        state=molecule("de novo primed template", (
            role("exposed genomic template"), role("7–9-nt DNA primers"))),
    )
    graph.react("prime", "Synthesize short DNA primers directly from dNTPs",
                StepKind.EXTENSION,
                ((template, InputRole.TEMPLATE), (primpol, InputRole.CATALYST)),
                (primed,), source=authority)
    phi29 = graph.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product(
        "amplified", "branched multiple-displacement WGA product",
        MaterialKind.MOLECULAR,
        state=molecule("TruePrime-amplified genomic DNA", (
            role("long genomic copies"), role("displaced single strands"))),
    )
    graph.react("amplify", f"Extend primers by strand displacement at 30 °C for {duration}",
                StepKind.EXTENSION,
                ((primed, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)),
                (amplified,), source=authority)
    heat = graph.input("inactivation_condition", "65 °C for 10 min", MaterialKind.CONDITION)
    product = graph.product(
        "inactivated_wga", "TruePrime-amplified genomic DNA after phi29 inactivation",
        MaterialKind.MOLECULAR, state=amplified.state)
    graph.react("inactivate", "Inactivate phi29 at 65 °C for 10 min",
                StepKind.THERMAL_TREATMENT,
                ((amplified, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (product,), source=authority)
    return product


def add_reported_sequencing(graph, amplified, *, authority: str):
    """Show the reported downstream library without assigning unknown NEBNext oligos."""
    condition = graph.input("fragmentation_condition", "downstream fragmentation condition",
                            MaterialKind.CONDITION, inferred=True, source=authority)
    fragments = graph.product("fragments", "fragmented TruePrime-amplified DNA",
                              MaterialKind.MOLECULAR, state=amplified.state,
                              inferred=True, source=authority)
    graph.react("fragment", "Fragment the amplified DNA for library preparation",
                StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (condition, InputRole.CONDITION)),
                (fragments,), source=authority + " (Covaris; library details unspecified)")
    library = graph.product("sequencing_library", "unspecified NEBNext library",
                            MaterialKind.LIBRARY, inferred=True, source=authority)
    graph.react("library_prep", "Prepare the unspecified NEBNext sequencing library",
                StepKind.MAGIC, ((fragments, InputRole.TEMPLATE),), (library,),
                reason="The paper names NEBNext but not its kit generation, reaction steps or adapter-end geometry.",
                source=authority)
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Read 1", "Read 2"),
        limitation="Paired-end HiSeq 2500 was reported; the kit version and primer sites are unknown."))
    return library
