"""Protocol-owned reaction graph for cropseq-multi."""
from types import SimpleNamespace
from multimodal_protocols import CELL, UMI, GUIDE, SAMPLE, _finish_rt_pcr
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from pcr import PCRGeometry

TITLE = "CROPseq-multi — multiplexed perturbation cassette"
DOI = "10.1101/2024.03.17.585235"
CITATION = "Walton et al., bioRxiv (2024)"
SUMMARY = ("A reverse-oriented tRNA–sgRNA array in the lentiviral 3′ LTR expresses multiple "
           "guides while short internal barcodes support mRNA, direct-guide and optical readouts.")
FINAL_CAPTION = ("A 12-base internal guide barcode remains linked to its spacer and can be "
                 "read from poly(A) mRNA, Pol III guide RNA or the CSMv2 T7-IVT product.")
NOTES = "01_cropseq_multi.html"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a>.'
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles",
                             (CELL, UMI, GUIDE, role("linked guide spacer")), index=SAMPLE)
SEQ_PRIMERS = ()
SEQUENCING_INTRO = FINAL_CAPTION
READ_LENGTHS = None
SEQUENCING_ENDING = FINAL_CAPTION
SEQUENCING_UNAVAILABLE = ("Target panels, pool barcodes and kit revisions vary; the defining "
                          "source establishes the read roles, while this page does not invent "
                          "a universal primer sequence.")
oligos = lambda: ()
_definition = SimpleNamespace(title=TITLE, doi=DOI, citation=CITATION, endpoint=FINAL_CAPTION)

def _cropseq_multi():
    return chain(molecule("CROPseq-multi 3-prime LTR cassette", (
        role("U6 promoter, reverse orientation"), role("tRNA"), role("sgRNA 1"),
        role("tRNA"), role("sgRNA 2"), GUIDE)), (
        ("Duplicate the engineered 3′ LTR during lentiviral integration",
         molecule("integrated perturbation cassette", (role("5-prime LTR copy"),
             role("selection-gene transcript"), role("3-prime LTR copy"))),
         "Keeping the Pol III cassette on the minus strand prevents tRNA processing from cleaving lentiviral genomic RNA or the Pol II selection transcript."),
        ("Transcribe and tRNA-process the multiplexed Pol III guide array",
         molecule("two mature perturbation guides", (role("sgRNA 1 with iBAR"),
                                                       role("sgRNA 2 with iBAR"))),
         "RNases P and Z release separate guides; each 12-base iBAR lies only 19 bp from its linked spacer."),
        ("Capture perturbation identity by mRNA, direct guide capture or T7-IVT",
         molecule("CROPseq-multi readout", (CELL, UMI, GUIDE,
                                             role("linked guide spacer"))),
         "The same iBAR is transcribed antisense in polyadenylated mRNA and within the Pol III sgRNA; CSMv2 additionally includes a T7 promoter."),
    ))

workflow = _cropseq_multi


def _cropseq_multi_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cassette = graph.input("cassette", "CROPseq-multi lentiviral 3-prime-LTR cassette",
                           MaterialKind.MOLECULAR, state=_cropseq_multi().initial, primary=True)
    cells = graph.input("target_cells", "target cells", MaterialKind.BIOLOGICAL)
    integrase = graph.input("integrase", "lentiviral integration machinery",
                            MaterialKind.CATALYST)
    integrated = graph.product("integrated_cassette", "integrated duplicated-LTR perturbation cassette",
                               MaterialKind.MOLECULAR,
                               state=molecule("integrated perturbation cassette", (
                                   role("5-prime LTR copy"), role("selection-gene transcript"),
                                   role("3-prime LTR copy"))))
    graph.react("integrate", "Integrate and duplicate the engineered 3-prime LTR",
                StepKind.INTEGRATION,
                ((cassette, InputRole.TEMPLATE), (cells, InputRole.ANALYTE),
                 (integrase, InputRole.CATALYST)), (integrated,))
    rna_polymerase = graph.input("rna_polymerase", "Pol III transcription machinery",
                                 MaterialKind.CATALYST)
    guides = graph.product("guides", "tRNA-processed guides carrying linked iBARs",
                           MaterialKind.MOLECULAR,
                           state=molecule("processed perturbation guides", (
                               role("sgRNA 1 with iBAR"), role("sgRNA 2 with iBAR"))))
    graph.react("transcribe_guides", "Transcribe and tRNA-process the multiplexed guide array",
                StepKind.TRANSCRIPTION,
                ((integrated, InputRole.TEMPLATE), (rna_polymerase, InputRole.CATALYST)),
                (guides,))
    mrna = graph.product("mrna_branch", "poly(A)-mRNA perturbation-identity branch",
                         MaterialKind.MOLECULAR, state=guides.state)
    direct = graph.product("direct_branch", "direct-guide-capture branch",
                           MaterialKind.MOLECULAR, state=guides.state)
    t7 = graph.product("t7_branch", "CSMv2 T7-IVT branch", MaterialKind.MOLECULAR,
                       state=guides.state)
    graph.react("split_readouts", "Route perturbation identity to the three supported readouts",
                StepKind.SPLIT, ((guides, InputRole.TEMPLATE),), (mrna, direct, t7),
                changes_molecule=False)
    reviewed = ["cassette", "target_cells", "integrase", "integrate", "rna_polymerase",
                "transcribe_guides", "split_readouts"]
    mrna_final = role_library("CROPseq-multi mRNA-readout roles",
                              (CELL, UMI, GUIDE, role("linked guide spacer")), index=SAMPLE)
    direct_final = role_library("CROPseq-multi direct-guide roles",
                                (CELL, UMI, GUIDE, role("linked guide spacer")), index=SAMPLE)
    _finish_rt_pcr(graph, "mrna", mrna, mrna_final, "mRNA branch: " + d.endpoint,
                   reviewed, label="CROPseq-multi mRNA")
    _finish_rt_pcr(graph, "direct", direct, direct_final,
                   "direct-guide branch: " + d.endpoint, reviewed,
                   label="CROPseq-multi direct guide")
    t7_polymerase = graph.input("t7_rna_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    t7_rna = graph.product("t7_rna", "T7-amplified CSMv2 guide RNA", MaterialKind.MOLECULAR,
                           state=t7.state)
    graph.react("t7_transcription", "Amplify the CSMv2 guide-barcode product by T7 IVT",
                StepKind.TRANSCRIPTION,
                ((t7, InputRole.TEMPLATE), (t7_polymerase, InputRole.CATALYST)), (t7_rna,))
    reviewed.extend(("t7_rna_polymerase", "t7_transcription"))
    _finish_rt_pcr(graph, "t7", t7_rna, final, "T7-IVT branch: " + d.endpoint,
                   reviewed, label="CROPseq-multi T7-IVT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    """Assemble the source-scoped material graph here in this protocol directory."""
    return _cropseq_multi_graph(_definition, FINAL_LIBRARY)
