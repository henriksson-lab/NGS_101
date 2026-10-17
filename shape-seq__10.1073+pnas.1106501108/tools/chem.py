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

# SI p.1 prints the 58-nt A_adapter_t; P5 and the Read 1 site share ACAC, present once.
# (Model correction 2026-10-17: the earlier il.P5 + il.TRUSEQ_READ1 duplicated it, 62 nt.)
SHAPE_A_T = il.P5 + il.TRUSEQ_READ1[4:]

SHAPE_A_B = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
SHAPE_PE_READ2 = il.PE_READ2_PRIMER

# SI p.1: RT primers carry a 5' handle RRRY for (+) and YYYR for (-) reactions; an inline
# degenerate condition tag read first in Read 1, not a designed whitelist barcode.
SHAPE_HANDLE_FEATURE = feature("shape_treatment_handle", "inline_barcode", "degenerate",
                               note="RRRY = 1M7-treated (+); YYYR = solvent control (-)")

SHAPE_RNA_FEATURE = feature("shape_rna_identity", "feature_barcode", "whitelist",
                            whitelist="designed RNA identity barcodes")

SHAPE_FINAL_LIBRARY = Construct((
    Segment("P5", il.P5, "p5"),
    Segment("A_adapter_t / Read 1 site remainder", SHAPE_A_T[len(il.P5):], "r1"),
    Segment("plus/minus handle", "RRRY", "cbc", placeholder=True,
            feature=SHAPE_HANDLE_FEATURE),
    # SI Fig. S1 / Table S1: the 3' structure cassette ends 5'-linker hairpin-[4-nt bar
    # code]-RT primer binding site-3'.  The RT primer therefore extends straight into
    # the bar code, so the bar code sits next to the RT primer (Read 1 end), and the
    # SHAPE stop lies at the opposite, A_adapter_b-ligated end (main text: the two
    # pieces of information "are on opposite ends of the cDNA molecules").
    role("target-specific RT primer (cassette RT-primer binding site)"),
    role("cDNA of RNA identity barcode", tag="cbc", length_bp=4,
         molecular_feature=SHAPE_RNA_FEATURE),
    role("cDNA of 3-prime linker hairpin"),
    role("cDNA ending before the SHAPE-modified nucleotide"),
    Segment("A_adapter_b sequencing arm", SHAPE_A_B[:-len(il.P7_RC)], "r2"),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="SHAPE-Seq v1 library")


TITLE = 'SHAPE-Seq — multiplexed RNA flexibility probing'
NOTES = '01_shape_seq.html'
SOURCE = 'Lucks et al., PNAS (2011), <a href="https://doi.org/10.1073/pnas.1106501108">doi:10.1073/pnas.1106501108</a>.'
SUMMARY = 'SHAPE chemistry acylates flexible riboses, reverse transcriptase stops at the adducts, and paired-end sequencing links each stop to a pre-attached RNA identity barcode.'
CAVEAT = 'The original SHAPE-Seq protocol assays designed RNA pools in vitro; later in-cell and transcriptome-wide SHAPE methods use different capture architectures.'
FINAL_LIBRARY = SHAPE_FINAL_LIBRARY
FINAL_CAPTION = 'Source-exact v1 adapter arms flank a target-specific RT product. Read 1 begins with the four-base treatment handle, then the RT primer and the cDNA copy of the 3′-cassette RNA bar code; Read 2 starts at the A_adapter_b ligation junction and locates the SHAPE-dependent RT stop. Model correction (2026-10-17): an earlier version drew the RNA bar code beside A_adapter_b; SI Fig. S1/Table S1 place it between the 3′ linker and the RT-primer binding site. The target-specific RT-primer bases are not printed in the SHAPE-Seq SI.'
SEQ_PRIMERS = (
    sp.TRUSEQ['R1'],
    sp.custom('Read 2', 'Illumina PE Read 2', SHAPE_PE_READ2,
              'Illumina 11 January 2011 adapter and sequencing-primer letter, p. 6'),
)
READ_LENGTHS = {'Read 1': 50, 'Read 2': 50}
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining method amplified the library for 9 or 12 PCR cycles (SI; the main text reports nine) with primers binding the Illumina adapter sequences (bases not printed), without size selection, and used 2×50-cycle paired-end sequencing on a Genome Analyzer IIx.'
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return (f'A_adapter_t / PE_forward — 5′-{SHAPE_A_T}-3′', f'A_adapter_b — 5′-{SHAPE_A_B}/3SpC3/-3′ (5′ phosphorylation is inferred for ligation, not printed in the SHAPE-Seq SI)', 'treated RT primer pattern — 5′-[A_adapter_t]-RRRY-[target-specific RT primer]-3′', 'control RT primer pattern — 5′-[A_adapter_t]-YYYR-[target-specific RT primer]-3′')

def workflow() -> Workflow:
    return chain(molecule("bar-coded in-vitro transcribed RNA pool", (role("folded RNA"),
                  role("3-prime linker hairpin"), role("4-nt RNA bar code"),
                  role("RT-primer binding site"))), (
        ("Modify flexible riboses with 1M7 (+) or treat with solvent only (-)",
         molecule("SHAPE-modified RNA", (role("RNA with 2-prime-O-adduct at flexible nucleotide"),
                  role("3-prime structure cassette with bar code"))),
         "SI: 6.5 mM 1M7; the bar code was transcribed from a DNA template made with degenerate primers, so many RNAs share the chemistry and library workflow."),
        ("Primer-extend with RRRY (+) or YYYR (-) handle RT primers until reverse transcriptase stops",
         molecule("SHAPE-stop cDNA", (role("A_adapter_t"), role("treatment handle"),
                  role("RT primer"), role("cDNA of RNA bar code"),
                  role("cDNA ending before modified nucleotide"))),
         "(+) and (-) pools are kept separate; treated and control stop distributions are interpreted jointly."),
        ("Hydrolyse the RNA with NaOH (95 °C, 5 min)",
         molecule("single-stranded SHAPE-stop cDNA", (role("A_adapter_t"), role("treatment handle"),
                  role("RT primer"), role("cDNA of RNA bar code"), role("stop-defined cDNA 3-prime end"))),
         "SI p.1: RNA hydrolysed with 4 M NaOH, cDNA ethanol precipitated."),
        ("Ligate 3-prime-C3-blocked A_adapter_b to the cDNA 3-prime end with CircLigase",
         molecule("dual-adapter SHAPE cDNA", (role("A_adapter_t"), role("treatment handle"),
                  role("RT primer"), role("cDNA of RNA bar code"),
                  role("stop-defined cDNA"), role("A_adapter_b"))),
         "SI p.1: 68 °C 6 h, separate (+) and (-) ligations, 80 °C 10 min stop; excess adapter removed with AMPure XP."),
        ("PCR-amplify (9 or 12 cycles) and paired-end sequence",
         molecule("SHAPE-Seq library", (role("P5"), role("treatment handle"),
                  role("RNA bar code"), role("SHAPE-stop coordinate"), role("P7"))),
         "Maximum-likelihood comparison of plus- and minus-reagent channels yields nucleotide reactivities."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('modify', 'Modify flexible riboses with 1M7 (or solvent control)', StepKind.CHEMICAL_CONVERSION), ('rt', 'Primer-extend to SHAPE adducts with handle-tailed RT primers', StepKind.REVERSE_TRANSCRIPTION), ('hydrolyse', 'Hydrolyse the RNA template with NaOH', StepKind.CHEMICAL_CONVERSION), ('ligate', 'Ligate 3-prime-C3-blocked A_adapter_b with CircLigase', StepKind.LIGATION))

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
