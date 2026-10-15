"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('scTEM-seq', '10.1038/s41598-022-09765-x', 'Hunt et al., Scientific Reports (2022)', 'Bisulfite-converted single-cell DNA is amplified with indexed SINE Alu or LINE-1 primers while RNA from the same cell is processed in parallel.', 'Paired reads span targeted transposable-element CpGs; internal and adapter indexes identify the cell library.', '01_sctem_seq.html')
PAYLOAD = (CELL, role("internal dual index"), role("TE amplicon"))

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
    return chain(molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))), (
        ("Separate genomic DNA and mRNA for parallel processing",
         molecule("matched DNA and RNA branches", (CELL, role("genomic DNA"), role("full-length cDNA"))),
         "The split permits Smart-seq-type transcriptomics alongside targeted methylation."),
        ("Bisulfite-convert genomic DNA",
         molecule("converted DNA", (role("converted SINE Alu / LINE-1 templates"),)),
         "Unmethylated cytosines convert while methylated CpGs remain cytosine."),
        ("Amplify transposable elements with staggered, internally indexed primer pairs",
         molecule("scTEM amplicons", (CELL, role("internal dual index"), role("SINE Alu insert"), SAMPLE)),
         "A 0–5 N stagger and reversed primer orientation diversify the low-complexity amplicon pool."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one sorted cell containing genomic DNA and mRNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))),
        primary=True)
    dna = graph.product("dna_branch", "single-cell genomic DNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    rna = graph.product("rna_branch", "matched single-cell mRNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_modalities", "Separate genomic DNA and mRNA",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_dna", "bisulfite-converted transposable-element templates",
        MaterialKind.MOLECULAR,
        state=molecule("converted DNA", (
            role("converted SINE Alu / LINE-1 templates"),)))
    graph.react("bisulfite_conversion", "Bisulfite-convert genomic DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    te_forward = graph.input("te_forward", "staggered internally indexed TE forward-primer pool",
                             MaterialKind.OLIGO)
    te_reverse = graph.input("te_reverse", "staggered internally indexed TE reverse-primer pool",
                             MaterialKind.OLIGO)
    te_polymerase = graph.input("te_polymerase", "transposable-element PCR polymerase",
                                MaterialKind.CATALYST)
    dna_library = graph.product("dna_library", "scTEM transposable-element amplicon library",
                                MaterialKind.LIBRARY, state=final)
    graph.react("te_pcr", "Amplify converted transposable elements",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (te_forward, InputRole.STRUCTURAL),
                 (te_reverse, InputRole.STRUCTURAL), (te_polymerase, InputRole.CATALYST)),
                (dna_library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The panel uses staggered internally indexed primer families rather than one universal amplicon."))
    graph.sequence(dna_library, SequencingPlan(limitation="DNA branch: " + d.endpoint))
    rna_final = role_library(
        "scTEM matched transcriptome library roles",
        (CELL, role("full-length transcript cDNA")), index=SAMPLE)
    reviewed = [
        "cell", "separate_modalities", "bisulfite", "bisulfite_conversion",
        "te_forward", "te_reverse", "te_polymerase", "te_pcr", "dna_library",
    ]
    _rt_then_pcr(graph, "rna", rna, rna_final,
                 "matched RNA branch: Smart-seq-type full-length transcript reads",
                 reviewed, rt_primer_label="Smart-seq-type oligo(dT) RT primer",
                 label="scTEM transcriptome")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
