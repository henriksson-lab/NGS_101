"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('MulTI-Tag', '10.1038/s41587-022-01522-9', 'Meers et al., Nature Biotechnology (2022)', 'Sequential antibody-targeted pA–Tn5 reactions insert target-specific adapter barcodes before single-cell combinatorial indexing.', 'Each genomic fragment carries both protein-target identity and cell identity.', '01_multi_tag.html')
PAYLOAD = (CELL, TARGET, role("target-proximal genomic insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ('Published Supplementary Table 1 prints the target-barcoded pA–Tn5 '
          'adapter and amplification-primer families; the source manifest confirms '
          'that publisher MOESM2 and preprint media-2 are byte-identical. This '
          'role-level diagram has not yet assembled a specific target/index '
          'combination into a final molecule.')
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("permeabilized cell with chromatin targets", (role("target 1"), role("target 2"), role("target 3"))), (
        ("Bind the first primary antibody and its barcoded pA–Tn5 reagent",
         molecule("first target tagged", (TARGET, role("target-1 chromatin fragment"))),
         "The transposon barcode identifies the antibody target."),
        ("Repeat antibody-directed tagmentation sequentially for additional targets",
         molecule("multiplex tagged chromatin", (TARGET, role("target-proximal genomic DNA"))),
         "Sequential binding and stringent washes retain target-specific adapter identities in one cell."),
        ("Combinatorially index single cells and amplify",
         molecule("MulTI-Tag library", (CELL, TARGET, role("target-proximal genomic insert"), SAMPLE)),
         "Target barcode and cell barcode are independent structured identifiers."),
    ))


def oligos():
    return ('Tn5MErev — 5′-/5Phos/CTGTCTCTTATACACATCT-3′',
            'P5_i5_1 (H3K27me3 example) — 5′-TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG-3′')


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cells", "permeabilized cells containing several chromatin targets",
        MaterialKind.BIOLOGICAL,
        state=molecule("permeabilized cell", (
            role("target 1"), role("target 2"), role("target 3"))), primary=True)
    reviewed = ["cells"]

    def targeting_round(prefix: str, template, label: str):
        antibody = graph.input(prefix + "_antibody", label + " primary antibody",
                               MaterialKind.SELECTION)
        bound = graph.product(prefix + "_bound", label + " antibody-bound chromatin",
                              MaterialKind.MOLECULAR, state=template.state)
        unbound = graph.product(prefix + "_unbound", label + " unbound chromatin",
                                MaterialKind.WASTE)
        graph.react(prefix + "_bind", "Bind " + label + " primary antibody",
                    StepKind.BINDING,
                    ((template, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                    (bound, unbound))
        graph.discard(unbound, "does not retain the selected antibody-target complex")
        transposome = graph.input(prefix + "_transposome", label + " barcode-loaded pA–Tn5",
                                  MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", label + " pA–Tn5 transposase",
                          MaterialKind.CATALYST)
        tagged = graph.product(prefix + "_tagged", label + " target-barcoded chromatin",
                               MaterialKind.MOLECULAR,
                               state=molecule(label + " tagged chromatin", (
                                   TARGET, role("target-proximal genomic DNA"))))
        losses = graph.product(prefix + "_losses", label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(prefix + "_tagment", "Tagment " + label + " target-proximal DNA",
                    StepKind.TAGMENTATION,
                    ((bound, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                     (tn5, InputRole.CATALYST)), (tagged, losses))
        graph.discard(losses, "does not retain the target-specific adapter identity")
        reviewed.extend((prefix + "_antibody", prefix + "_bind",
                         prefix + "_transposome", prefix + "_tn5", prefix + "_tagment"))
        return tagged

    first = targeting_round("first", cells, "first target")
    multiplex = targeting_round("additional", first, "additional-target pool")
    partition = graph.input("cell_index_partition", "single-cell combinatorial-index partitioning condition",
                            MaterialKind.CONDITION)
    indexed = graph.product("indexed_cells", "combinatorially indexed target-tagged cells",
                            MaterialKind.MOLECULAR, state=multiplex.state)
    collisions = graph.product("index_collisions", "ambiguous cell-index partitions",
                               MaterialKind.WASTE)
    graph.react("index_cells", "Combinatorially index single cells",
                StepKind.PARTITION,
                ((multiplex, InputRole.TEMPLATE), (partition, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not yield an unambiguous cell-barcode combination")
    reviewed.extend(("cell_index_partition", "index_cells"))
    _finish_pcr(graph, "library", indexed, final, d.endpoint, reviewed,
                label="MulTI-Tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
