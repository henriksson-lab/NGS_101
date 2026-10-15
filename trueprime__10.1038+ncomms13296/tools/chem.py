"""Protocol-owned source, oligos, workflow and reaction graph."""
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
        ("Lyse the cell, denature DNA and neutralize", template,
         "The single-cell workflow lyses on ice before neutralization; purified DNA is alkaline-denatured."),
        ("TthPrimPol synthesizes short DNA primers directly from dNTPs", primed,
         "The principal primer products are 7–9 nt; no synthetic random primer is added."),
        ("phi29 DNA polymerase extends the primers and displaces downstream strands", product,
         "Displacement exposes new single-stranded template for another priming round."),
        ("Repeat priming and strand-displacement synthesis at 30 °C", product,
         "The paper used 3 h for WGA and 3 or 6 h for single cells, then 10 min at 65 °C."),
    ))


FINAL_LIBRARY = Construct([role("TruePrime-amplified genomic DNA")], name="TruePrime WGA product")
SEQ_PRIMERS = ()
TITLE = 'TruePrime — primase-driven whole-genome amplification'
NOTES = '01_trueprime.html'
SOURCE = '<a href="https://doi.org/10.1038/ncomms13296">Picher et al., <i>Nature Communications</i> (2016)</a>.'
SUMMARY = 'TthPrimPol makes 7–9-nt DNA primers directly on denatured genomic DNA; phi29 DNA polymerase extends them by multiple-displacement synthesis.'
CAVEAT = 'TruePrime ends at amplified DNA. The paper used an unspecified NEBNext library preparation downstream, so no adapter sequence or sequencing-primer landing site is asserted here.'
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
    primpol = graph.input("tth_primpol", "TthPrimPol", MaterialKind.CATALYST)
    primed = graph.product("primed", "templates bearing de-novo 7–9-nt DNA primers",
                           MaterialKind.MOLECULAR, state=single.state)
    graph.react("prime", "Synthesize short DNA primers directly from dNTPs",
                StepKind.EXTENSION,
                ((single, InputRole.TEMPLATE), (primpol, InputRole.CATALYST)), (primed,))
    phi29 = graph.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "branched multiple-displacement WGA product",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("amplify", "Extend primers with repeated strand displacement",
                StepKind.EXTENSION,
                ((primed, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (amplified,))
    fragmentation = graph.input("fragmentation_condition", "downstream fragmentation condition",
                                MaterialKind.CONDITION, inferred=True,
                                source=authority + " (downstream NEBNext implementation unspecified)")
    fragments = graph.product("fragments", "fragmented TruePrime-amplified DNA",
                              MaterialKind.MOLECULAR, state=amplified.state, inferred=True,
                              source=authority + " (downstream NEBNext implementation unspecified)")
    graph.react("fragment", "Fragment the amplified DNA for library preparation",
                StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragments,), source=authority + " (downstream NEBNext kit unspecified)")
    adapter = graph.input("nebnext_adapter", "unspecified NEBNext sequencing adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=authority + " (kit generation and oligos not identified)")
    ligase = graph.input("adapter_ligase", "NEBNext adapter ligase",
                         MaterialKind.CATALYST, inferred=True,
                         source=authority + " (kit generation not identified)")
    library = graph.product("sequencing_library", "unspecified NEBNext library",
                            MaterialKind.LIBRARY, inferred=True,
                            source=authority + " (adapter architecture unavailable)")
    graph.react("library_prep", "Prepare the unspecified NEBNext sequencing library",
                StepKind.LIGATION,
                ((fragments, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper does not identify the NEBNext kit generation or adapter-end geometry."),
                source=authority + " (downstream kit named only)")
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Read 1", "Read 2"),
        limitation="The paper reports paired-end HiSeq 2500 sequencing but not the NEBNext kit generation, adapter oligos or primer sites."))
    graph.declare_reviewed_scope([
        "genomic_dna", "denaturation_condition", "denature", "tth_primpol", "prime",
        "phi29", "amplify", "fragmentation_condition", "fragment", "nebnext_adapter",
        "adapter_ligase", "library_prep", "sequencing_library"], authority)
    return graph.require_publishable()
