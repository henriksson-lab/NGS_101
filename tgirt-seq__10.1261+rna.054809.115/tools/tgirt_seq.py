"""Source-owned Qin2016 TGIRT-seq ligation architecture (original Fig. S1)."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'lib')]
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

TITLE = 'TGIRT-seq — Qin 2016 RNA-adapter / cDNA-ligation method'
NOTES = '01_tgirt-seq.html'
SOURCE = '<a href="https://doi.org/10.1261/rna.054809.115">Qin et al., <i>RNA</i> (2016)</a>, Methods, Figure 1 and original Supplemental Figure S1.'
SUMMARY = 'A blocked R2 RNA/R2R DNA starter initiates TGIRT template switching. RNA hydrolysis leaves single-stranded cDNA; a preadenylated R1R adapter is ligated to its 3′ end before single-index PCR.'
CAVEAT = ('This page represents the original 2016 method. Its starter and PCR oligos differ '
          'from the later 2021 UMI/dual-index protocol. The source prints ordered R1R '
          'with 5′ phosphate but describes the ligating input as 5′ adenylated; preparation '
          'of that activated input is an explicit unresolved step in this model. '
          'No selected barcode, barcode width or loaded instrument-primer reagent is assigned.')

# Source assignment is local; the shared Illumina constants supply only canonical bases.
R2_RNA = revcomp(il.TRUSEQ_READ2).replace('T', 'U')
R2_RNA_WRITTEN = ''.join('r' + base for base in R2_RNA) + '/3SpC3/'
R2R_DNA = il.TRUSEQ_READ2 + 'N'
R1R_DNA = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAG'
R1R_WRITTEN = '/5Phos/' + R1R_DNA + '/3SpC3/'
R1_SITE = revcomp(R1R_DNA)
R2_SITE = revcomp(il.TRUSEQ_READ2)
P5_PRIMER = il.P5 + il.TRUSEQ_SMALL_RNA_RA5
if il.P5[:23] + R1_SITE != P5_PRIMER:
    raise ValueError('Original R1R does not supply the source P5-primer annealing arm')
INDEX_FEATURE = feature('sample_index_i7', 'sample_index', 'unknown')
INDEX = role('sample index; width unassigned', tag='cbc', molecular_feature=INDEX_FEATURE)
INDEX_RC = role('sample-index reverse complement; width unassigned', tag='cbc', molecular_feature=INDEX_FEATURE)
P7_PARTS = (Segment('P7 PCR tail', il.P7, 'p7'), INDEX,
            Segment('P7-primer annealing arm', il.TRUSEQ_READ2, 'r2'))

# The overhang is the first base of the cDNA, paired to the target RNA's 3′ base.
# It must not be discarded or treated as a UMI/extra adapter base.
CDNA_PARTS = (Segment('R2R complementary core', il.TRUSEQ_READ2, 'r2'),
              Segment('initiating target-complementary base', 'N'),
              role('copied target cDNA beyond initiating base'))
ADAPTER_CDNA = Construct(CDNA_PARTS, name='single-stranded R2R-linked cDNA')
TWO_ADAPTER_CDNA = Construct((*CDNA_PARTS, Segment('3-prime R1R adapter', R1R_DNA, 'r1')),
                             name='5-prime R2R / cDNA / R1R 3-prime')
FINAL_LIBRARY = Construct((Segment('P5 PCR tail', il.P5[:23], 'p5'),
    Segment('source R1R reverse complement', R1_SITE, 'r1'),
    role('RNA-derived insert, including starter-paired terminal base'),
    Segment('R2R reverse complement', R2_SITE, 'r2'), INDEX_RC,
    Segment('P7 reverse complement', il.P7_RC, 'p7')), name='Qin2016 single-index library')
FINAL_CAPTION = 'PCR-derived strand orientation: P5 — R1R reverse complement — RNA-derived insert — R2R reverse complement — sample-index reverse complement — P7 reverse complement. Selected index and width remain unassigned.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_ENDING = 'The original plasma libraries were sequenced with paired 100-nt HiSeq2500 reads or paired 75-/150-nt NextSeq500 reads. Table1/FigureS4 distinguish datasets; no single run is silently chosen.'
SEQUENCING_UNAVAILABLE = 'The loaded sequencing-primer oligos and selected index are not declared in this model. PCR-binding sites are known; this limitation is not a claim that all linked sources omit instrument-primer bases.'
READ_SPECS = (
    sp.read_spec('Read 1', '100 HiSeq; 75 or 150 NextSeq', 'RNA-derived insert, source paired-end plasma libraries; no UMI is assigned to this original branch', 'include'),
    sp.read_spec('Index 1 (i7)', None, 'sample index; selected sequence and width unassigned', 'demultiplex'),
    sp.read_spec('Index 2 (i5)', 'not assigned', 'original S1 P5 has no i5 sample-index window; acquired technical cycles are not assigned', 'undetermined'),
    sp.read_spec('Read 2', '100 HiSeq; 75 or 150 NextSeq', 'opposite end of RNA-derived insert, source paired-end plasma libraries', 'include'),
)

def oligos():
    return ('R2 RNA (Fig. S1; 3′ spacer) — 5′-' + R2_RNA_WRITTEN + '-3′',
            'R2R DNA (equimolar single-base overhang) — 5′-' + R2R_DNA + '-3′',
            'R1R as ordered in Fig. S1 — 5′-' + R1R_WRITTEN + '-3′',
            'Illumina multiplex PCR primer — 5′-' + P5_PRIMER + '-3′',
            'Illumina barcode PCR primer — 5′-' + il.P7 + '[BARCODE]' + il.TRUSEQ_READ2 + '-3′')

def starter_scene(*, extended=False):
    if len(R2R_DNA) != len(R2_RNA) + 1 or revcomp(R2R_DNA[:-1]) != R2_RNA.replace('U', 'T'):
        raise ValueError('Original R2/R2R must form a full duplex with exactly one DNA 3-prime overhang')
    sc = Scene()
    sc.strand('rna', [Segment('R2 RNA', R2_RNA, 'rna')], label='blocked R2 RNA', mod3='/3SpC3/')
    parts = CDNA_PARTS if extended else CDNA_PARTS[:2]
    sc.anneal('dna', list(parts), to='rna', pair=('R2R complementary core', 'R2 RNA'),
              label='extended R2R DNA' if extended else 'R2R DNA primer')
    sc.note('dna', 'N is the target-pairing initiating base; not a UMI')
    return sc

def adapter_scene(*, activated):
    sc = Scene()
    sc.strand('adapter', [Segment('R1R adapter', R1R_DNA, 'r1')],
              label='ligation-ready R1R' if activated else 'ordered R1R',
              mod5='/5App/' if activated else '/5Phos/', mod3='/3SpC3/')
    if activated:
        sc.ligation_end('adapter', 'R1R adapter', kind='single-strand', polarity='5-prime', owner='adapter')
    return sc

def ligated_scene():
    sc = Scene()
    sc.strand('cdna', list(TWO_ADAPTER_CDNA), label='5′ R2R — cDNA — R1R 3′', mod3='/3SpC3/')
    sc.junction('cdna', CDNA_PARTS[-1].name, '3-prime R1R adapter', 'sealed cDNA 3′–R1R 5′ junction')
    return sc

def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    target = graph.input('target_rna', 'plasma RNA with accessible 3-prime hydroxyl', MaterialKind.BIOLOGICAL,
                          state=molecule('target RNA', (role('target RNA'),)), primary=True)
    starter = graph.input('starter_duplex', 'original R2 RNA / R2R DNA starter', MaterialKind.OLIGO,
                           state=starter_scene())
    tgirt = graph.input('tgirt', 'source dataset-specific TeI4c or GsI-IIC/TGIRT-III RT', MaterialKind.CATALYST)
    extended = graph.product('extended_dna', 'extended R2R DNA with RNA still present', MaterialKind.MOLECULAR,
                              state=starter_scene(extended=True))
    graph.react('template_switch', 'Template-switch and extend at 60°C for 15 min', StepKind.TEMPLATE_SWITCHING,
                ((target, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL), (tgirt, InputRole.CATALYST)), (extended,))
    alkali = graph.input('alkali', 'NaOH hydrolysis followed by HCl neutralization', MaterialKind.CONDITION)
    cdna = graph.product('completed_cdna', 'purified single-stranded R2R-linked cDNA', MaterialKind.MOLECULAR, state=ADAPTER_CDNA)
    waste = graph.product('hydrolyzed_rna', 'hydrolyzed RNA', MaterialKind.WASTE,
                          state=molecule('RNA hydrolysis products', (role('RNA fragments'),)))
    graph.react('remove_rna', 'Hydrolyze RNA and purify cDNA', StepKind.CHEMICAL_CONVERSION,
                ((extended, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)), (cdna, waste))
    graph.discard(waste, 'RNA is hydrolyzed before cDNA purification and ligation.')
    ordered = graph.input('ordered_r1r', 'Fig. S1 ordered 5-prime-phosphorylated R1R', MaterialKind.OLIGO,
                           state=adapter_scene(activated=False))
    adapter = graph.product('activated_r1r', 'Methods ligation-ready 5-prime-adenylated R1R', MaterialKind.ADAPTER,
                            state=adapter_scene(activated=True))
    graph.react('activate_r1r', 'Prepare the source-described preadenylated R1R', StepKind.MAGIC,
                ((ordered, InputRole.TEMPLATE),), (adapter,),
                reason='Fig. S1 prints 5Phos R1R; Methods requires 5-prime adenylated R1R. The activation procedure is not declared in this model; no particular adenylation enzyme or reagent is inferred.')
    ligase = graph.input('app_ligase', 'Thermostable 5-prime AppDNA/RNA Ligase (NEB)', MaterialKind.CATALYST)
    dual = graph.product('two_adapter_cdna', 'R2R–cDNA–R1R ligation product', MaterialKind.MOLECULAR, state=ligated_scene())
    graph.react('ligate_r1r', 'Ligate R1R 5-prime end to cDNA 3-prime hydroxyl', StepKind.LIGATION,
                ((cdna, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (dual,),
                geometry=AdapterLigationGeometry.known('single-strand', adapter_scene(activated=True),
                    statement='Original Methods and Fig.1 specify a 5-prime-adenylated, 3-prime-blocked R1R donor ligated to the cDNA 3-prime end.'))
    p5 = graph.input('p5_primer', 'original unindexed P5 PCR primer', MaterialKind.OLIGO,
                      state=Construct((Segment('P5 PCR primer', P5_PRIMER),), name='original P5 primer'))
    p7 = graph.input('p7_primer', 'original indexed P7 PCR primer; barcode unassigned', MaterialKind.OLIGO,
                      state=molecule('P7 indexed primer', P7_PARTS))
    polymerase = graph.input('pcr_polymerase', 'Phusion High-Fidelity DNA polymerase', MaterialKind.CATALYST)
    library = graph.product('library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react('index_pcr', 'Amplify with original P5/indexed P7 primers (12 Phusion cycles)', StepKind.PCR,
                ((dual, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(TWO_ADAPTER_CDNA,
                    PCRPrimer('P7 primer known 3-prime annealing arm; index tail unassigned', il.TRUSEQ_READ2, anneal_nt=len(il.TRUSEQ_READ2)),
                    PCRPrimer('source P5 PCR primer', P5_PRIMER, anneal_nt=len(R1_SITE)),
                    statement='Both source-known annealing arms face inward on the R2R–cDNA–R1R template. The indexed P7 primer’s 5-prime barcode tail is not assigned.'))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(tuple(graph.materials) + tuple(graph.steps), SOURCE)
    return graph.require_publishable()
