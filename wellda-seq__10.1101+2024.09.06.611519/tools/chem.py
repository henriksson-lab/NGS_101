"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('wellDA-seq', '10.1101/2024.09.06.611519', 'Wang et al., bioRxiv (2024)', 'Indexed accessibility tagging and a whole-genome branch co-profile chromatin accessibility and copy-number variation from the same single cells.', 'Matched ATAC and genomic reads share a single-cell well identity.', '01_wellda_seq.html')
PAYLOAD = (CELL, role("ATAC or genomic insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("single nucleus in a well", (role("genomic DNA"), role("accessible chromatin"))), (
        ("Tagment accessible chromatin with indexed Tn5",
         molecule("ATAC-marked genome", (CELL, role("accessible-DNA fragments"), role("remaining genome"))),
         "Tagmentation records accessibility and supplies the first molecular handle."),
        ("Displace the tagged products and amplify the residual genome",
         molecule("paired products from one nucleus", (CELL, role("ATAC branch"), role("whole-genome branch"))),
         "The branches preserve their common well identity while measuring accessibility and copy number."),
        ("Add branch-compatible Illumina indexes and pool wells",
         molecule("wellDA-seq libraries", (CELL, SAMPLE, role("ATAC or genomic insert"))),
         "The same single-cell well barcode links the epigenome and genome."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nucleus = graph.input(
        "nucleus", "single nucleus in a well",
        MaterialKind.BIOLOGICAL,
        state=molecule("single nucleus in a well", (
            role("genomic DNA"), role("accessible chromatin"))), primary=True)
    indexed_adapter = graph.input("indexed_transposome", "well-indexed loaded Tn5 adapter",
                                  MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    marked = graph.product(
        "marked_genome", "ATAC-marked genome retaining residual genomic DNA",
        MaterialKind.MOLECULAR,
        state=molecule("ATAC-marked genome", (
            CELL, role("accessible-DNA fragments"), role("remaining genome"))))
    losses = graph.product("tn5_losses", "non-library Tn5 products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin with indexed Tn5",
                StepKind.TAGMENTATION,
                ((nucleus, InputRole.TEMPLATE), (indexed_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (marked, losses))
    graph.discard(losses, "does not retain an amplifiable indexed accessibility fragment")
    atac = graph.product("atac_branch", "displaced ATAC-tagged products",
                         MaterialKind.MOLECULAR, state=marked.state)
    genome = graph.product("genome_branch", "residual genomic DNA",
                           MaterialKind.MOLECULAR, state=marked.state)
    graph.react("separate_branches", "Separate tagged products from the residual genome",
                StepKind.SPLIT, ((marked, InputRole.TEMPLATE),), (atac, genome),
                changes_molecule=False)
    reviewed = ["nucleus", "indexed_transposome", "tn5", "tagment",
                "separate_branches"]
    _finish_pcr(graph, "atac", atac, final, "ATAC branch: " + d.endpoint,
                reviewed, label="wellDA ATAC branch")
    random_primer = graph.input("genome_random_primer", "whole-genome amplification primer",
                                MaterialKind.OLIGO)
    wga_polymerase = graph.input("genome_wga_polymerase", "strand-displacing DNA polymerase",
                                 MaterialKind.CATALYST)
    amplified = graph.product("amplified_genome", "whole-genome-amplified residual genome",
                              MaterialKind.MOLECULAR, state=genome.state)
    graph.react("whole_genome_amplification", "Amplify the residual genome",
                StepKind.EXTENSION,
                ((genome, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (wga_polymerase, InputRole.CATALYST)), (amplified,))
    reviewed.extend(("genome_random_primer", "genome_wga_polymerase",
                     "whole_genome_amplification"))
    _finish_pcr(graph, "genome", amplified, final, "whole-genome branch: " + d.endpoint,
                reviewed, label="wellDA genomic branch")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
