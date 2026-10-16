"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
import seqprimers as sp
from rna_end_readout import (StopReadout, BaseConversionReadout, SelectiveCleavageReadout,
                             _finish_pcr, _apply_linear_step)

OTTR_PRIMER = Segment("universal +1T DNA primer", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT", "r2")

OTTR_DUPLEX_RNA = Segment("blocked RNA complement", "GATCGGAAGAGCACACGTCTGAACTCCAGTC", "r2")

OTTR_ACCEPTOR = Segment("universal adaptor template", "ACACTCTTTCCCTACACGACGCTCTTCCGATCTC", "r1")


TITLE = 'OTTR-seq — Ordered Two-Template Relay'
NOTES = '01_ottr_seq.html'
SOURCE = 'Upton et al., PNAS (2021); <a href="https://doi.org/10.1073/pnas.2107900118">doi:10.1073/pnas.2107900118</a>.'
SUMMARY = 'BoMoC standardizes RNA 3′ ends with one ddR, then performs two ordered template jumps to copy the insert and the second adapter in one reverse-transcription reaction.'
CAVEAT = 'The oligos displayed here are printed in the authors’ bioRxiv v1 Table S1 (PDF text lines 768–1200). The final PNAS article cites its SI Table S1, but that appendix is currently inaccessible; exact agreement between versions remains unresolved.'
FINAL_LIBRARY = Construct((role('P5'), role('i5', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'random')), role('R1 adapter'), role('input cDNA'), role('R2 adapter'), role('i7', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'random')), role('P7')), name='OTTR-seq library')
FINAL_CAPTION = 'Universal-adapter branch after indexed PCR; the paper also validates a full-length, PCR-free branch.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The paper sequenced OTTR libraries on an Illumina MiniSeq 75-cycle kit.'
SEQUENCING_UNAVAILABLE = 'The paper specifies Illumina adapter families but does not publish a distinct custom sequencing primer.'

def oligos():
    return ('universal +1T DNA primer — 5′-/5Cy5/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT-3′', 'blocked RNA complement — 5′-rGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrU/3SpC3/-3′', 'universal adapter template — 5′-/5Cy5/ACACTCTTTCCCTACACGACGCTCTTCCrGrArUrCrUrC-3′')

def workflow() -> Workflow:
    return chain(molecule("input small RNA", (role("RNA of any end sequence"),)), (
        ("Cap each input 3′ end with one ddA and/or ddG",
         molecule("3′-blocked input", (role("input RNA"), role("single nonextendable ddR"))),
         "BoMoC terminal-transferase activity runs in Mn²⁺; phosphatase then destroys free ddRTP."),
        ("Switch BoMoC from the +1Y primer duplex onto input RNA",
         molecule("first template-jump product", (OTTR_PRIMER, role("input-derived cDNA"))),
         "A +1T/+1C DNA end pairs with the standardized ddA/ddG RNA end."),
        ("Add one dG and relay onto the 3′C adapter template",
         molecule("dual-adapter cDNA", (OTTR_PRIMER, role("input-derived cDNA"), OTTR_ACCEPTOR)),
         "The second, ordered template jump copies one terminating adapter template."),
        ("PCR-add flow-cell ends and indexes",
         molecule("OTTR-seq library", (role("P5 and i5"), role("R1 adapter"), role("input cDNA"),
                  role("R2 adapter"), role("i7 and P7"))),
         "Universal-adapter libraries use four to eight PCR cycles; full-length adapters can omit PCR."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('cap_end', 'Cap each input 3-prime end with one ddA and/or ddG', StepKind.TAILING), ('first_jump', 'Template-switch from the standardized primer duplex onto input RNA', StepKind.TEMPLATE_SWITCHING), ('second_jump', 'Relay onto the terminating adapter template', StepKind.TEMPLATE_SWITCHING))

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in RECIPE:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed)
    _finish_pcr(
        graph, "library", current, FINAL_LIBRARY,
        SEQUENCING_ENDING or FINAL_CAPTION, reviewed, label=TITLE,
        geometry=PCR_GEOMETRY,
        primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS,
        read_uses=READ_USES)
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
