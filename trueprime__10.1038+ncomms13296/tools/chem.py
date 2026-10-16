"""Protocol-owned source, oligos, workflow and reaction graph."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import ThreePrimeOverhangAdapter, ligation_product_scene, react_dA_tail_then_ligate, repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow
from protocol_families.primpol_mda import (AMPLIFICATION_C, INACTIVATION_C,
                                          PRIMER_LENGTH_NT, add_amplification,
                                          add_reported_sequencing)

TRUEPRIME_PRIMER_LENGTH_NT = PRIMER_LENGTH_NT

TRUEPRIME_AMPLIFICATION_C = AMPLIFICATION_C

TRUEPRIME_INACTIVATION_C = INACTIVATION_C


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
CAVEAT = 'This page follows the purified-DNA WGA experiment. The <a href="../trueprime-single-cell__10.1038+ncomms13296/trueprime-single-cell.html">single-cell L2-lysis workflow</a> has a separate page. TruePrime ends at amplified DNA. The paper used an unspecified NEBNext library preparation downstream, so no adapter sequence or sequencing-primer landing site is asserted here.'
FINAL_CAPTION = 'The named protocol product is amplified genomic DNA, before sequencing adapters.'
SEQUENCING_INTRO = ''
READ_LENGTHS = None
SEQUENCING_ENDING = 'The paper fragmented this product, prepared a NEBNext library and performed paired-end HiSeq 2500 sequencing.'
SEQUENCING_UNAVAILABLE = 'The publication does not identify the NEBNext kit generation or publish its adapter and sequencing-primer sequences.'
def oligos():
    return ()

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    authority = SOURCE
    graph = ProtocolGraph("TruePrime — primase-driven whole-genome amplification",
                          authority=authority)
    dna = graph.input("genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    denaturation = graph.input("denaturation_condition", "alkaline denaturation and neutralization",
                               MaterialKind.CONDITION)
    single = graph.product("single_stranded", "denatured single-stranded genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("denatured genomic DNA", (role("single-stranded genomic template"),)))
    graph.react("denature", "Denature genomic DNA and neutralize",
                StepKind.THERMAL_TREATMENT,
                ((dna, InputRole.TEMPLATE), (denaturation, InputRole.CONDITION)), (single,))
    amplified = add_amplification(graph, single, authority=authority, duration="3 h")
    add_reported_sequencing(graph, amplified, authority=authority)
    graph.declare_reviewed_scope([
        "genomic_dna", "denaturation_condition", "denature", "tth_primpol", "prime",
        "phi29", "amplify", "inactivation_condition", "inactivate", "fragmentation_condition", "fragment",
        "library_prep", "sequencing_library"], authority)
    return graph.require_publishable()
