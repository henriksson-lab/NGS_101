"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('Butt-seq', '10.1101/gad.350434.123', 'Yu and Rosbash, Genes & Development (2023)', 'TGIRT template switching captures native nascent-RNA 3′ termini without RNA fragmentation, followed by second-strand synthesis and indexed PCR.', 'The terminal insert base reports the nascent RNA 3′ end; a separate UMI read supports molecule deduplication.', '01_butt_seq.html')
PAYLOAD = (UMI, role("nascent-RNA terminal insert"))

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
    return chain(molecule("purified chromatin-associated nascent RNA", (role("nascent RNA 5-prime body"), role("native RNA 3-prime end"))), (
        ("TGIRT template-switch from an RNA/DNA starter duplex onto the RNA 3′ end",
         molecule("first-strand cDNA", (UMI, role("nascent RNA 3-prime-end cDNA"))),
         "The first copied base reports the polymerase-proximal RNA end without RNA fragmentation."),
        ("Switch salt concentration and synthesize the second strand",
         molecule("double-stranded Butt-seq product", (UMI, role("nascent-RNA insert"))),
         "The salt switch favors the intended first template switch and suppresses later spurious switching."),
        ("PCR-amplify, size-select and index",
         molecule("Butt-seq library", (UMI, role("nascent-RNA terminal insert"), SAMPLE)),
         "A separate UMI read supports molecule-aware paired-end analysis."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "nascent_rna", "purified chromatin-associated nascent RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("purified chromatin-associated nascent RNA", (
            role("nascent RNA 5-prime body"), role("native RNA 3-prime end"))),
        primary=True)
    starter = graph.input("starter_duplex", "UMI-bearing TGIRT RNA/DNA starter duplex",
                          MaterialKind.OLIGO)
    tgirt = graph.input("tgirt", "TGIRT reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product(
        "first_strand", "nascent-RNA 3-prime-end first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA", (
            UMI, role("nascent RNA 3-prime-end cDNA"))))
    graph.react("template_switch", "Template-switch TGIRT onto the native RNA 3-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    dna_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                 MaterialKind.CATALYST)
    double = graph.product(
        "double_stranded", "double-stranded Butt-seq product",
        MaterialKind.MOLECULAR,
        state=molecule("double-stranded Butt-seq product", (
            UMI, role("nascent-RNA insert"))))
    graph.react("second_strand", "Switch salt and synthesize the second strand",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (dna_polymerase, InputRole.CATALYST)), (double,))
    reviewed = ["nascent_rna", "starter_duplex", "tgirt", "template_switch",
                "second_strand_primer", "second_strand_polymerase", "second_strand"]
    forward = graph.input("pcr_forward", "Butt-seq PCR forward primer", MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "indexed Butt-seq PCR reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified Butt-seq library",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "PCR-amplify and index", StepKind.PCR,
                ((double, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines the UMI and read roles without one universal complete primer geometry."))
    size_rule = graph.input("size_selection", "Butt-seq library-size selection rule",
                            MaterialKind.SELECTION)
    library = graph.product("sequencing_library", "size-selected Butt-seq library",
                            MaterialKind.LIBRARY, state=final)
    rejected = graph.product("rejected_sizes", "off-size Butt-seq products",
                             MaterialKind.WASTE)
    graph.react("select_size", "Size-select the indexed library", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (library, rejected))
    graph.discard(rejected, "falls outside the selected nascent-RNA library size")
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
                     "size_selection", "select_size", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
