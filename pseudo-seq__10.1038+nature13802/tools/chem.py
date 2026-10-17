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

PSEUDO_ADAPTER = Segment("3-prime adapter; 5-prime chemistry unresolved", "TGGAATTCTCGGGTGCCAAGG", "r2")

PSEUDO_RT = Segment("Pseudo-seq RT primer", "GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATTCACTCAGCCTTGGCACCCGAGAATTCCA", "r2")

PSEUDO_FWD = Segment("PCR forward", "AATGATACGGCGACCACCGA", "p5")

PSEUDO_REV = Segment("indexed PCR reverse", "CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", "p7", placeholder=True,
                     feature=feature("sample_index", "sample_index", "unknown"))

# The printed reverse primer split into its tail, barcode window and RT-primer-matching end.
PSEUDO_REV_STATE = Construct((
    Segment("P7 tail", PSEUDO_REV.top[:24], "p7"),
    Segment("XXXXXX sample barcode", PSEUDO_REV.top[24:30], "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown")),
    Segment("RT-primer-matching 3-prime end", PSEUDO_REV.top[30:], "r2"),
), name="barcoded PCR reverse primer")


TITLE = 'Pseudo-seq — transcriptome-wide pseudouridine mapping'
NOTES = '01_pseudo_seq.html'
SOURCE = 'Carlile et al., Nature (2014); <a href="https://doi.org/10.1038/nature13802">doi:10.1038/nature13802</a>.'
SUMMARY = 'Persistent Ψ–CMC adducts stop reverse transcriptase; enriched cDNA termini relative to a mock-treated control locate pseudouridines.'
CAVEAT = ("<b>Published adapter discrepancy.</b> Carlile et al., Methods, ‘Pseudo-seq Library Preparation’ (PMC4224642, paragraph on 3′ ligation) call the adapter pre-adenylated, but print its IDT sequence as <code>/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/</code>. A 5′ phosphate and a 5′ adenylate are different states. The same paragraph runs the ligation with T4 RNA ligase (NEB M0204) in 1X buffer without ATP, which is consistent with an adenylated donor but does not settle the ordered modification. This page follows the printed base sequence and modification string while depicting ligation as described; whether the ordered or reacted adapter was adenylated remains unresolved.")
FINAL_LIBRARY = Construct((role('P5'), role('stop-defined RNA cDNA'), role('sample index', tag='cbc', length_bp=6, molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='Pseudo-seq library')
FINAL_CAPTION = 'PCR product from circularized, stop-defined cDNA.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining experiment sequenced the amplified library on an Illumina HiSeq 2000.'
SEQUENCING_UNAVAILABLE = 'The paper publishes library oligos but does not identify a separate custom sequencing-primer oligo.'

def oligos():
    return ('3′ adapter — 5′-/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/-3′', 'RT primer — 5′-/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA-3′', 'PCR forward — 5′-AATGATACGGCGACCACCGA-3′', 'indexed PCR reverse — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′')

def workflow() -> Workflow:
    readout = StopReadout("CMC", "pseudouridine", "parallel mock −CMC")
    return chain(molecule("fragmented RNA", (role("RNA containing U and Ψ"),)), (  # ZnCl2 (yeast) or Zn acetate (HeLa) fragmentation precedes CMC
        ("Derivatize RNA with CMC and reverse U/G adducts",
         molecule("CMC-treated RNA", (role("RNA"), role("persistent Ψ–CMC adduct"))),
         f"{readout.control} distinguishes Ψ-dependent stops from background."),
        ("Repair ends, size-select, and ligate the reported 3′ adapter",
         molecule("adapter-ligated RNA", (role("CMC-treated RNA fragment"), PSEUDO_ADAPTER),
                  junctions=(("CMC-treated RNA fragment", "3-prime adapter; 5-prime chemistry unresolved", "RNA ligation"),)),
         "CIP/PNK treatment precedes 8% urea-PAGE selection of 80–100, 100–120 and 120–140-nt RNA fragments; T4 RNA ligase 1 joins the adapter in ATP-free buffer."),
        ("Reverse-transcribe until AMV RT stops one nucleotide 3′ of Ψ–CMC",
         molecule("truncated cDNA", (PSEUDO_RT, role("cDNA ending at Ψ-dependent stop"))),
         "Truncated cDNAs are gel size-selected; the 5′ end of each read reports the position of the block."),
        ("Circularize cDNA and amplify with indexed primers",
         molecule("Pseudo-seq library", (PSEUDO_FWD, role("stop-defined cDNA"), PSEUDO_REV)),
         "CircLigase closes the single-stranded cDNA before Illumina PCR."),
    ))


PCR_GEOMETRY = None
PCR_PRIMER_SOURCE = ("Carlile et al. (2014), Methods, Pseudo-seq Library Preparation "
                     "(PMC4224642): forward primer and barcoded reverse primer.")
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('fragment', 'Fragment RNA with zinc at high temperature', StepKind.FRAGMENTATION), ('cmc', 'Derivatize pseudouridine with persistent CMC', StepKind.CHEMICAL_CONVERSION), ('repair', 'Dephosphorylate RNA fragment ends with CIP and PNK', StepKind.DEPHOSPHORYLATION), ('ligate', 'Ligate the reported 3-prime adapter; 5-prime chemistry unresolved', StepKind.LIGATION), ('rt', 'Reverse-transcribe to the pseudouridine-dependent stop', StepKind.REVERSE_TRANSCRIPTION), ('circle', 'Circularize the stop-defined cDNA', StepKind.CIRCULARIZATION))

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
        step_label="Amplify circularized cDNA with the barcoded reverse primer (Phusion)",
        forward_label="PCR forward primer",
        reverse_label="barcoded PCR reverse primer (XXXXXX sample barcode)",
        polymerase_label="Phusion DNA polymerase (NEB M0530)",
        forward_state=PSEUDO_FWD, reverse_state=PSEUDO_REV_STATE,
        primer_source=PCR_PRIMER_SOURCE,
        geometry=PCR_GEOMETRY,
        primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS,
        read_uses=READ_USES)
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
