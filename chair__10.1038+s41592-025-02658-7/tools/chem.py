"""Protocol-owned reaction graph for chair."""
from types import SimpleNamespace
from multimodal_protocols import CELL, UMI, SAMPLE, _finish_pcr, _finish_rt_pcr
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from pcr import PCRGeometry

TITLE = "ChAIR — chromatin accessibility, interaction and RNA"
DOI = "10.1038/s41592-025-02658-7"
CITATION = "Chai et al., Nature Methods (2025)"
SUMMARY = ("In-situ restriction digestion and proximity ligation record chromatin contacts, "
           "Tn5 marks accessible anchors, and 10x Multiome barcodes the DNA and RNA branches together.")
FINAL_CAPTION = ("ChAIR-RNA, ChAIR-ATAC and ChAIR-PET share a droplet cell barcode; paired "
                 "DNA ends distinguish contact junctions from accessibility fragments.")
NOTES = "01_chair.html"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a>.'
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles",
                             (CELL, UMI, role("RNA, ATAC or contact insert")), index=SAMPLE)
SEQ_PRIMERS = ()
SEQUENCING_INTRO = FINAL_CAPTION
READ_LENGTHS = None
SEQUENCING_ENDING = FINAL_CAPTION
SEQUENCING_UNAVAILABLE = ("Target panels, pool barcodes and kit revisions vary; the defining "
                          "source establishes the read roles, while this page does not invent "
                          "a universal primer sequence.")
oligos = lambda: ()
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=FINAL_CAPTION)

def _chair():
    return chain(molecule("crosslinked nucleus", (role("chromatin contacts"),
                                                   role("accessible chromatin"),
                                                   role("poly(A) RNA"))), (
        ("Digest chromatin in situ with a restriction enzyme",
         molecule("restricted crosslinked chromatin", (role("restriction fragment A"),
                                                        role("crosslinked restriction fragment B"))),
         "Restriction cleavage exposes ends while crosslinks preserve spatial neighbours."),
        ("Proximity-ligate crosslinked restriction ends",
         molecule("3D contact junction", (role("restriction fragment A"),
                                           role("restriction fragment B")),
                  junctions=(("restriction fragment A", "restriction fragment B", "proximity ligation"),)),
         "The ligation junction records a chromatin contact."),
        ("Tagment accessible chromatin while retaining contact products",
         molecule("ChAIR DNA precursor", (role("Tn5 end"), role("accessible anchor"),
                                           role("contact partner"), role("Tn5 end"))),
         "The DNA branch yields ordinary accessibility fragments and paired-end contact tags anchored at accessible chromatin."),
        ("Bar-code DNA and RNA together in 10x Multiome droplets",
         molecule("three matched readouts", (CELL, UMI, role("ChAIR-RNA"),
                                              role("ChAIR-ATAC"), role("ChAIR-PET"))),
         "One droplet cell barcode links gene expression, accessibility and proximity-ligation junctions."),
    ))

workflow = _chair


def _chair_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nucleus = graph.input("nucleus", "crosslinked nucleus containing chromatin and poly(A) RNA",
                          MaterialKind.BIOLOGICAL, state=_chair().initial, primary=True)
    restriction = graph.input("restriction_enzyme", "in-situ restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("digested", "restricted crosslinked chromatin",
                             MaterialKind.MOLECULAR,
                             state=molecule("restricted chromatin", (
                                 role("restriction fragment A"),
                                 role("crosslinked restriction fragment B"), role("poly(A) RNA"))))
    graph.react("digest", "Digest chromatin in situ", StepKind.DIGESTION,
                ((nucleus, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "proximity DNA ligase", MaterialKind.CATALYST)
    contacts = graph.product("contacts", "chromatin-contact junctions retained with RNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("3D contact junction", (
                                 role("restriction fragment A"), role("restriction fragment B"),
                                 role("poly(A) RNA")),
                                 junctions=(("restriction fragment A", "restriction fragment B",
                                             "proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate crosslinked restriction ends",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    transposome = graph.input("transposome", "ATAC loaded transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "accessible-anchor and contact products retained with RNA",
                           MaterialKind.MOLECULAR, state=_chair().reactions[2].after)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin while retaining contact products",
                StepKind.TAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain an amplifiable accessible or contact fragment")
    droplets = graph.input("droplets", "10x Multiome droplet-barcoding condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell-barcoded RNA, ATAC and contact products",
                             MaterialKind.MOLECULAR,
                             state=molecule("three matched readouts", (CELL, UMI,
                                 role("ChAIR-RNA"), role("ChAIR-ATAC"), role("ChAIR-PET"))))
    failed = graph.product("failed_partitions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Barcode DNA and RNA together in Multiome droplets",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable cell barcode")
    rna = graph.product("rna_branch", "ChAIR-RNA branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    atac = graph.product("atac_branch", "ChAIR-ATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    pet = graph.product("pet_branch", "ChAIR-PET contact branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    graph.react("split_readouts", "Separate RNA, accessibility and contact readouts",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (rna, atac, pet),
                changes_molecule=False)
    reviewed = ["nucleus", "restriction_enzyme", "digest", "proximity_ligase",
                "proximity_ligation", "transposome", "tn5", "tagment", "droplets",
                "partition", "split_readouts"]
    rna_final = role_library("ChAIR-RNA roles", (CELL, UMI, role("transcript cDNA")), index=SAMPLE)
    atac_final = role_library("ChAIR-ATAC roles", (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    pet_final = role_library("ChAIR-PET roles", (CELL, role("contact-junction insert")), index=SAMPLE)
    _finish_rt_pcr(graph, "rna", rna, rna_final, "ChAIR-RNA: " + d.endpoint,
                   reviewed, label="ChAIR RNA")
    _finish_pcr(graph, "atac", atac, atac_final, "ChAIR-ATAC: " + d.endpoint,
                reviewed, label="ChAIR ATAC")
    _finish_pcr(graph, "pet", pet, pet_final, "ChAIR-PET: " + d.endpoint,
                reviewed, label="ChAIR contact")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    """Assemble the source-scoped material graph here in this protocol directory."""
    return _chair_graph(_definition, FINAL_LIBRARY)
