"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('SPEAC-seq', '10.1038/s41596-024-01056-1', 'Faust Akl et al., Nature Protocols (2025)', 'CRISPR-perturbed and reporter cells are paired and cultured in droplets; functional fluorescence sorting precedes guide-cassette PCR and sequencing.', 'Guide counts are compared between functionally sorted and control droplet fractions; this is a pooled screen, not a transcriptome library.', '01_speac_seq.html')
PAYLOAD = (GUIDE, role("guide-flanking amplicon"))

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
    return chain(molecule("pooled CRISPR-perturbed effector cells", (GUIDE, role("reporter cells"))), (
        ("Pair one perturbed cell with one reporter cell in a picolitre droplet",
         molecule("co-culture droplet", (GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))),
         "Different dyes gate one-to-one pairing before the interaction develops."),
        ("Sort droplets by the functional fluorescent response",
         molecule("phenotype-enriched droplet fraction", (GUIDE, role("reporter state"))),
         "The screen's phenotype is a sorted population, not a single-cell transcriptome."),
        ("Break droplets, recover genomic DNA and PCR-amplify the integrated guide cassette",
         molecule("SPEAC-seq library", (GUIDE, SAMPLE, role("guide-flanking amplicon"))),
         "Guide counts in response-positive and control fractions identify interaction regulators."),
    ))


def oligos():
    return ("lentiCRISPRv2 P7 index primer — 5′-CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT-3′",)


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "cell_pools", "CRISPR-perturbed effector cells and fluorescent reporter cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("pooled screen cells", (GUIDE, role("effector cells"),
                                                role("reporter cells"))), primary=True)
    droplets = graph.input("pairing_condition", "picolitre one-to-one droplet pairing condition",
                           MaterialKind.CONDITION)
    paired = graph.product("paired_droplets", "one-effector/one-reporter droplets",
                           MaterialKind.MOLECULAR,
                           state=molecule("co-culture droplet", (
                               GUIDE, role("perturbed cell"), role("fluorescent reporter cell"))))
    pairing_failures = graph.product("pairing_failures", "empty or incorrectly paired droplets",
                                     MaterialKind.WASTE)
    graph.react("pair_cells", "Pair one perturbed and one reporter cell per droplet",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (paired, pairing_failures))
    graph.discard(pairing_failures, "does not contain the required one-to-one cell pair")
    response_gate = graph.input("response_gate", "functional fluorescence sort gates",
                                MaterialKind.SELECTION)
    positive = graph.product("positive_droplets", "response-positive droplet fraction",
                             MaterialKind.MOLECULAR, state=paired.state)
    control = graph.product("control_droplets", "response-control droplet fraction",
                            MaterialKind.MOLECULAR, state=paired.state)
    graph.react("sort_response", "Sort droplets by functional response",
                StepKind.SELECTION,
                ((paired, InputRole.TEMPLATE), (response_gate, InputRole.SELECTION)),
                (positive, control))
    reviewed = ["cell_pools", "pairing_condition", "pair_cells", "response_gate",
                "sort_response"]

    def recover(prefix: str, template, label: str):
        extraction = graph.input(prefix + "_extraction", label + " droplet-break and DNA-extraction condition",
                                 MaterialKind.CONDITION)
        dna = graph.product(prefix + "_dna", label + " integrated guide-cassette DNA",
                            MaterialKind.MOLECULAR,
                            state=molecule(label + " guide DNA", (
                                GUIDE, role("guide-flanking genomic DNA"))))
        graph.react(prefix + "_recover", "Break droplets and recover genomic DNA",
                    StepKind.EXTRACTION,
                    ((template, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                    (dna,))
        reviewed.extend((prefix + "_extraction", prefix + "_recover"))
        _finish_pcr(graph, prefix, dna, final, label + ": " + d.endpoint,
                    reviewed, label=label + " guide-cassette")

    recover("positive", positive, "response-positive")
    recover("control", control, "control")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
