"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('LiMCA', '10.1038/s41592-024-02239-0', 'Wu et al., Nature Methods (2024)', 'Cytoplasm and nucleus are separated: Smart-seq2 reads full-length RNA while restriction digestion, proximity ligation and META amplification recover high-depth chromatin contacts.', 'Matched full-length transcript and Hi-C libraries are joined by the recorded cell of origin.', '01_limca.html')
PAYLOAD = (role("Smart-seq2 insert or Hi-C junction"),)

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
    return chain(molecule("one intact cell", (role("cytoplasmic mRNA"), role("nuclear chromatin"))), (
        ("Physically separate cytoplasm and nucleus",
         molecule("matched compartments", (CELL, role("cytoplasm"), role("nucleus"))),
         "Compartment separation lets each modality retain its high-sensitivity chemistry."),
        ("Amplify cytoplasmic RNA by Smart-seq2",
         molecule("full-length transcriptome branch", (CELL, role("full-length cDNA"))),
         "The RNA branch reports gene expression from the same cell."),
        ("Crosslink, restriction-digest and proximity-ligate the nucleus",
         molecule("single-cell contact junctions", (role("restriction fragment A"), role("restriction fragment B")),
                  junctions=(("restriction fragment A", "restriction fragment B", "proximity ligation"),)),
         "The nuclear branch uses conventional 3C chemistry."),
        ("Apply META whole-genome amplification and build the Hi-C library",
         molecule("paired LiMCA libraries", (CELL, role("Smart-seq2 insert or Hi-C junction"), SAMPLE)),
         "Transcript and contact libraries are integrated by their shared cell-of-origin record."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one intact cell",
        MaterialKind.BIOLOGICAL,
        state=molecule("one intact cell", (
            role("cytoplasmic mRNA"), role("nuclear chromatin"))), primary=True)
    cytoplasm = graph.product("cytoplasm", "cytoplasmic RNA compartment",
                              MaterialKind.MOLECULAR, state=cell.state)
    nucleus = graph.product("nucleus", "matched nuclear chromatin compartment",
                            MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_compartments", "Physically separate cytoplasm and nucleus",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (cytoplasm, nucleus),
                changes_molecule=False)
    reviewed = ["cell", "separate_compartments"]
    rna_final = role_library("LiMCA Smart-seq2 transcriptome library roles",
                             (CELL, role("full-length transcript cDNA")), index=SAMPLE)
    _rt_then_pcr(graph, "rna", cytoplasm, rna_final,
                 "RNA branch: " + d.endpoint, reviewed,
                 rt_primer_label="Smart-seq2 oligo(dT) primer",
                 label="LiMCA Smart-seq2")

    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_nucleus", "crosslinked nuclear chromatin",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    graph.react("crosslink", "Crosslink nuclear chromatin", StepKind.CROSS_LINKING,
                ((nucleus, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    restriction = graph.input("restriction_enzyme", "Hi-C restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("restriction_fragments", "restriction-digested fixed chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Restriction-digest the nucleus", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)),
                (digested,))
    ligase = graph.input("proximity_ligase", "in-nucleus proximity DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product(
        "contacts", "single-cell proximity-ligation junctions",
        MaterialKind.MOLECULAR,
        state=molecule("single-cell contact junctions", (
            role("restriction fragment A"), role("restriction fragment B")),
            junctions=(("restriction fragment A", "restriction fragment B",
                        "proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contacts,))
    reversal = graph.input("crosslink_reversal", "cross-link reversal condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified_contacts", "purified single-cell contact DNA",
                             MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover contact DNA",
                StepKind.DEPROTEINIZATION,
                ((contacts, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (purified,))
    meta_primer = graph.input("meta_primer", "META whole-genome amplification primer",
                              MaterialKind.OLIGO)
    meta_polymerase = graph.input("meta_polymerase", "META strand-displacing polymerase",
                                  MaterialKind.CATALYST)
    amplified = graph.product("amplified_contacts", "META-amplified contact DNA",
                              MaterialKind.MOLECULAR, state=purified.state)
    graph.react("meta_wga", "Apply META whole-genome amplification",
                StepKind.EXTENSION,
                ((purified, InputRole.TEMPLATE), (meta_primer, InputRole.STRUCTURAL),
                 (meta_polymerase, InputRole.CATALYST)), (amplified,))
    reviewed.extend((
        "formaldehyde", "crosslink", "restriction_enzyme", "digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "meta_primer", "meta_polymerase", "meta_wga",
    ))
    _finish_pcr(graph, "hic", amplified, final, "Hi-C branch: " + d.endpoint,
                reviewed, label="LiMCA Hi-C")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
