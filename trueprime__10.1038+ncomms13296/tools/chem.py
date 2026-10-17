"""Protocol-owned source, oligos, workflow and reaction graph."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from chemdraw import Construct, Workflow
from protocol_model import molecule, role
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow

# Picher 2016: principal primase products and purified-DNA WGA Methods.
TRUEPRIME_PRIMER_LENGTH_NT = (7, 9)
TRUEPRIME_AMPLIFICATION_C = 30
TRUEPRIME_INACTIVATION_C = 65


def workflow() -> Workflow:
    template = molecule("denatured genomic DNA", (role("single-stranded genomic template"),))
    primed = molecule("de novo primed template", (
        role("single-stranded genomic template"), role("7–9-nt DNA primer")))
    product = molecule("branched multiple-displacement product", (
        role("long genomic copies"), role("displaced single strands"), role("new DNA primers")))
    input_dna = molecule("genomic DNA", (role("genomic DNA"),))
    return _workflow(input_dna, (
        ("Alkaline-denature purified DNA for 3 min and neutralize", template,
         "The source uses Buffer D followed by Buffer N for the purified-DNA WGA experiment."),
        ("TthPrimPol synthesizes short DNA primers directly from dNTPs", primed,
         "The principal primer products are 7–9 nt; no synthetic random primer is added."),
        ("phi29 DNA polymerase extends the primers and displaces downstream strands", product,
         "Displacement exposes new single-stranded template for another priming round."),
        ("Repeat priming and strand-displacement synthesis for 3 h at 30 °C", product,
         "The purified-DNA WGA experiment ends with 10 min at 65 °C."),
    ))


FINAL_LIBRARY = Construct([role("TruePrime-amplified genomic DNA")], name="TruePrime WGA product")
SEQ_PRIMERS = ()
TITLE = 'TruePrime — purified-DNA whole-genome amplification'
NOTES = '01_trueprime.html'
SOURCE = '<a href="https://doi.org/10.1038/ncomms13296">Picher et al., <i>Nature Communications</i> (2016)</a>.'
SUMMARY = 'Purified genomic DNA is alkaline-denatured, then TthPrimPol makes 7–9-nt primers that phi29 DNA polymerase extends by strand-displacement synthesis.'
CAVEAT = ('This page follows the purified-DNA WGA experiment (Methods "Whole genome amplification": 3 h at 30 °C for Fig. 4b,d). '
          'The <a href="../trueprime-single-cell__10.1038+ncomms13296/trueprime-single-cell.html">single-cell L2-lysis workflow</a> has a separate page. '
          'TruePrime ends at amplified DNA. The only purified-DNA products the paper reports sequencing are the 1 ng, 1 pg and 1 fg inputs of '
          'Supplementary Fig. 2, amplified for 6 h and sequenced on Ion Torrent; their library preparation is not described. '
          'The Methods NEBNext/HiSeq 2500 passage and the ENA runs belong to the single-cell and non-amplified samples, so this page does not assign them to purified DNA.')
FINAL_CAPTION = 'The named protocol product is amplified genomic DNA, before sequencing adapters.'
SEQUENCING_INTRO = ''
READ_LENGTHS = None
SEQUENCING_ENDING = 'Supplementary Fig. 2 reports Ion Torrent reads from 6 h purified-DNA TruePrime products (1 ng, 1 pg, 1 fg inputs).'
SEQUENCING_UNAVAILABLE = 'The paper does not describe the Ion Torrent library preparation, adapters or read lengths for these products.'
def oligos():
    return ()

# Source-owned Picher 2016 amplification and downstream experiment.
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
    """Supplementary Fig. 2 Ion Torrent check; its library preparation is not described."""
    supp = authority + " Supplementary Fig. 2 legend"
    library = graph.product("sequencing_library", "Ion Torrent library with undescribed preparation",
                            MaterialKind.LIBRARY, inferred=True, source=supp)
    graph.react("library_prep", "Prepare an Ion Torrent library from 6 h TruePrime products",
                StepKind.MAGIC, ((amplified, InputRole.TEMPLATE),), (library,),
                reason=("Supplementary Fig. 2 states only that 1 ng, 1 pg and 1 fg purified-DNA inputs were "
                        "amplified for 6 h and sequenced on Ion Torrent; fragmentation, library kit and "
                        "adapters are not described."),
                source=supp)
    graph.sequence(library, SequencingPlan(
        limitation=("Ion Torrent reads used to show >95% human-assigned reads; read layout and length are not "
                    "reported. The 6 h amplification differs from the 3 h Fig. 4 Methods reaction drawn above.")))
    return library


def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    authority = SOURCE
    graph = ProtocolGraph("TruePrime — primase-driven whole-genome amplification",
                          authority=authority)
    dna = graph.input("genomic_dna", "2.5 µl purified human genomic DNA (6 pg, or 1 ng to 100 ag)", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    denaturation = graph.input("denaturation_condition", "alkaline denaturation and neutralization",
                               MaterialKind.CONDITION)
    single = graph.product("single_stranded", "denatured single-stranded genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("denatured genomic DNA", (role("single-stranded genomic template"),)))
    graph.react("denature", "Denature with 2.5 µl Buffer D for 3 min at room temperature, then neutralize with 2.5 µl Buffer N",
                StepKind.THERMAL_TREATMENT,
                ((dna, InputRole.TEMPLATE), (denaturation, InputRole.CONDITION)), (single,),
                source=authority + " Methods, Whole genome amplification")
    amplified = add_amplification(graph, single, authority=authority, duration="3 h")
    add_reported_sequencing(graph, amplified, authority=authority)
    graph.declare_reviewed_scope([
        "genomic_dna", "denaturation_condition", "denature", "tth_primpol", "prime",
        "phi29", "amplify", "inactivation_condition", "inactivate",
        "library_prep", "sequencing_library"], authority)
    return graph.require_publishable()
