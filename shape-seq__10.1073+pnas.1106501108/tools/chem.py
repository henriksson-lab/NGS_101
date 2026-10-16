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

SHAPE_A_T = il.P5 + il.TRUSEQ_READ1

SHAPE_A_B = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
SHAPE_PE_READ2 = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"

SHAPE_HANDLE_FEATURE = feature("shape_treatment_handle", "feature_barcode", "whitelist",
                               whitelist="RRRY treated / YYYR control handle families")

SHAPE_RNA_FEATURE = feature("shape_rna_identity", "feature_barcode", "whitelist",
                            whitelist="designed RNA identity barcodes")

SHAPE_FINAL_LIBRARY = Construct((
    Segment("P5", il.P5, "p5"),
    Segment("A_adapter_t / Read 1 site", il.TRUSEQ_READ1, "r1"),
    Segment("plus/minus handle", "RRRY", "cbc", placeholder=True,
            feature=SHAPE_HANDLE_FEATURE),
    role("target-specific RT-primer site"),
    role("SHAPE-stop-defined cDNA"),
    role("RNA identity barcode", tag="cbc", length_bp=4,
         molecular_feature=SHAPE_RNA_FEATURE),
    Segment("A_adapter_b sequencing arm", SHAPE_A_B[:-len(il.P7_RC)], "r2"),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="SHAPE-Seq v1 library")


TITLE = 'SHAPE-Seq — multiplexed RNA flexibility probing'
NOTES = '01_shape_seq.html'
SOURCE = 'Lucks et al., PNAS (2011), <a href="https://doi.org/10.1073/pnas.1106501108">doi:10.1073/pnas.1106501108</a>.'
SUMMARY = 'SHAPE chemistry acylates flexible riboses, reverse transcriptase stops at the adducts, and paired-end sequencing links each stop to a pre-attached RNA identity barcode.'
CAVEAT = 'The original SHAPE-Seq protocol assays designed RNA pools in vitro; later in-cell and transcriptome-wide SHAPE methods use different capture architectures.'
FINAL_LIBRARY = SHAPE_FINAL_LIBRARY
FINAL_CAPTION = 'Source-exact v1 adapter arms flank a target-specific RT product. The first paired read begins with the four-base treatment handle; the opposite read locates the SHAPE-dependent RT stop.'
SEQ_PRIMERS = (
    sp.TRUSEQ['R1'],
    sp.custom('Read 2', 'Illumina PE Read 2', SHAPE_PE_READ2,
              'Illumina 11 January 2011 adapter and sequencing-primer letter, p. 6'),
)
READ_LENGTHS = {'Read 1': 50, 'Read 2': 50}
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining method used 2×50-cycle paired-end sequencing on a Genome Analyzer IIx.'
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return (f'A_adapter_t / PE_forward — 5′-{SHAPE_A_T}-3′', f'A_adapter_b — 5′-{SHAPE_A_B}/3SpC3/-3′ (5′ phosphorylation is inferred for ligation, not printed in the SHAPE-Seq SI)', 'treated RT primer pattern — 5′-[A_adapter_t]-RRRY-[target-specific RT primer]-3′', 'control RT primer pattern — 5′-[A_adapter_t]-YYYR-[target-specific RT primer]-3′')

def workflow() -> Workflow:
    return chain(molecule("barcoded in-vitro RNA pool", (role("RNA-specific 3-prime DNA barcode"),
                  role("folded RNA"))), (
        ("Modify flexible riboses with a SHAPE reagent and retain a minus-reagent control",
         molecule("SHAPE-modified RNA", (role("RNA"),
                  role("2-prime-O-adduct at flexible nucleotide"))),
         "The pre-attached RNA barcode lets many defined RNAs share the chemical and library workflow."),
        ("Primer-extend until reverse transcriptase reaches a SHAPE adduct",
         molecule("SHAPE-stop cDNA", (role("common primer handle"),
                  role("cDNA ending before modified nucleotide"), role("RNA identity barcode"))),
         "Treated and control stop distributions are interpreted jointly rather than as presence/absence calls."),
        ("Ligate the second sequencing adapter to cDNA",
         molecule("dual-adapter SHAPE cDNA", (role("paired-end adapter"),
                  role("stop-defined cDNA"), role("RNA identity barcode"),
                  role("paired-end adapter"))),
         "The two ends separately identify the assayed RNA and the primer-extension stop."),
        ("PCR-amplify and paired-end sequence",
         molecule("SHAPE-Seq library", (role("P5"), role("SHAPE-stop coordinate"),
                  role("RNA identity barcode"), role("P7"))),
         "Maximum-likelihood comparison of plus- and minus-reagent channels yields nucleotide reactivities."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('modify', 'Modify flexible riboses with a SHAPE reagent', StepKind.CHEMICAL_CONVERSION), ('rt', 'Primer-extend to SHAPE adducts', StepKind.REVERSE_TRANSCRIPTION), ('ligate', 'Ligate the second sequencing adapter to cDNA', StepKind.LIGATION))

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
