"""Bounded paired-end library preparation for small chromatin footprints."""
from batch_ngs import truseq_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import da_tailed


def henikoff_single_index_library(name: str, payload):
    """Published short-insert paired-end shell shared by ChEC and ORGANIC."""
    return truseq_library(
        list(payload), name, dual_index=False, i7_nt=6)


def _finish_henikoff_library(graph: ProtocolGraph, prefix: str, template, label: str,
                              final, primers, source: str, reviewed: list[str],
                              read_lengths, adapter_geometry, p5_sequence: str, p7_sequence: str,
                              primer_source_statement: str) -> None:
    """Published low-input end repair through paired-end library construction."""
    repair = graph.input(prefix + "_repair_enzymes", label + " end-repair enzymes",
                         MaterialKind.CATALYST)
    repaired = graph.product(prefix + "_repaired", label + " end-repaired DNA",
                             MaterialKind.MOLECULAR, state=template.state)
    graph.react(prefix + "_repair", "End-repair " + label + " DNA",
                StepKind.END_REPAIR,
                ((template, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    kinase = graph.input(prefix + "_kinase", label + " T4 polynucleotide kinase",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product(prefix + "_phosphorylated", label + " phosphorylated DNA",
                                   MaterialKind.MOLECULAR, state=repaired.state)
    graph.react(prefix + "_phosphorylate", "Phosphorylate repaired " + label + " ends",
                StepKind.PHOSPHORYLATION,
                ((repaired, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    tailing = graph.input(prefix + "_tailing_enzyme", label + " Klenow exo-minus",
                          MaterialKind.CATALYST)
    tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                           MaterialKind.MOLECULAR,
                           state=da_tailed(label + " dA-tailed DNA", label + " insert"))
    graph.react(prefix + "_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((phosphorylated, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,))
    adapters = graph.input(prefix + "_adapters", label + " published indexed TruSeq adapters",
                           MaterialKind.ADAPTER)
    ligase = graph.input(prefix + "_adapter_ligase", label + " adapter ligase",
                         MaterialKind.CATALYST)
    adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final)
    graph.react(prefix + "_ligation", "Ligate indexed TruSeq adapters",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=adapter_geometry)
    p5 = graph.input(prefix + "_p5_primer", label + " published P5 enrichment primer",
                     MaterialKind.OLIGO)
    p7 = graph.input(prefix + "_p7_primer", label + " published P7 enrichment primer",
                     MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", label + " KAPA PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    geometry = PCRGeometry.exact(
        final, PCRPrimer("P5 enrichment primer", p5_sequence),
        PCRPrimer("P7 enrichment primer", p7_sequence),
        statement=primer_source_statement)
    graph.react(prefix + "_pcr", "PCR-enrich the " + label + " library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry)
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers), read_lengths))
    reviewed.extend((
        prefix + "_repair_enzymes", prefix + "_repair", prefix + "_kinase",
        prefix + "_phosphorylate", prefix + "_tailing_enzyme", prefix + "_tail",
        prefix + "_adapters", prefix + "_adapter_ligase", prefix + "_ligation",
        prefix + "_p5_primer", prefix + "_p7_primer", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))
