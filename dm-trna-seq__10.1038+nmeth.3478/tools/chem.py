"""Source-owned DM-TGIRT states; molecular primitives contain no protocol recipe."""
from __future__ import annotations

from chemdraw import Construct, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import closed_circle, molecule, role

TITLE = 'DM-tRNA-seq (DM-TGIRT-seq)'
DOI = '10.1038/nmeth.3478'
CITATION = 'Zheng et al., Nature Methods (2015)'
SUMMARY = 'Wild-type and D135S AlkB demethylation precedes TGIRT template switching, RNA hydrolysis, single-stranded cDNA circularization and indexed PCR.'
NOTES = '01_dm_trna_seq.html'
SOURCE = 'Zheng et al., Nature Methods (2015); <a href="https://doi.org/10.1038/nmeth.3478">doi:10.1038/nmeth.3478</a>.'
CAVEAT = ('INFERRED — The complementary DNA starter below is a minimal reconstruction: '
          'full pairing to the printed 41-nt RNA plus the source-specified single 3′ T. '
          'Additional DNA extensions and the position of the reported 32P label are not assigned. '
          'The final library remains a role-level endpoint; exact circular-template PCR assembly '
          'and the authors’ loaded HiSeq sequencing reagents are not declared in this model.')
SEQUENCING_ENDING = 'The 2015 Methods give both PCR primers and Illumina HiSeq sequencing. The loaded run-primer reagent is not assigned here.'
SEQUENCING_UNAVAILABLE = 'This model does not declare an exact circular-template PCR product or loaded HiSeq run-primer set; no exact sequencing landing geometry is asserted.'
SEQUENCING_INTRO = ''
SEQ_PRIMERS = ()
READ_LENGTHS = None
FINAL_CAPTION = 'Role-level library after the source-published PCR. Index width, selected barcode and complete product assembly are unassigned model data.'

RNA_SEQUENCE = 'AGAUCGGAAGAGCACACGUCUAGUUCUACAGUCCGACGAUC'
DM_TEMPLATE = Segment('blocked RNA starter', RNA_SEQUENCE, 'rna')
DM_FORWARD = Segment('multiplex PCR primer', il.P5 + il.TRUSEQ_SMALL_RNA_RA5, 'p5')
INDEX = role('sample index; width unassigned', tag='cbc',
             molecular_feature=feature('sample_index', 'sample_index', 'unknown'))
REVERSE_PARTS = (Segment('P7', il.P7, 'p7'), INDEX,
                 Segment('reverse PCR arm', il.TRUSEQ_READ2, 'r2'))
# RNA partitions into the two opposed PCR-primer annealing regions. The remaining
# primer bases are 5-prime tails, so full-primer containment is not required.
FORWARD_ANNEAL = il.TRUSEQ_SMALL_RNA_RA5[-20:]
REVERSE_ANNEAL = il.TRUSEQ_READ2[-21:]
if RNA_SEQUENCE.replace('U', 'T') != revcomp(REVERSE_ANNEAL) + FORWARD_ANNEAL:
    raise ValueError('DM starter no longer partitions into the published PCR binding regions')

def oligos():
    return ('TGIRT RNA template — 5′-' + RNA_SEQUENCE + '/3SpC3/-3′',
            'multiplex PCR primer — 5′-' + DM_FORWARD.top + '-3′',
            'indexed PCR primer — 5′-' + il.P7 + '[BARCODE]' + il.TRUSEQ_READ2 + '-3′')

def dna_parts():
    return (Segment('minimal complementary DNA core', revcomp(RNA_SEQUENCE), 'r1', inferred=True),
            Segment('source-specified initiating T', 'T'),
            role('tRNA-derived cDNA beyond initiating T'))

def starter_scene(*, extended=False):
    sc = Scene()
    sc.strand('rna', [DM_TEMPLATE], label='blocked RNA starter', mod3='/3SpC3/')
    parts = dna_parts() if extended else dna_parts()[:2]
    sc.anneal('dna', list(parts), to='rna',
              pair=('minimal complementary DNA core', DM_TEMPLATE.name),
              label='INFERRED — minimal DNA reconstruction')
    sc.note('dna', '32P label reported; position unassigned')
    return sc

FINAL_LIBRARY = Construct((role('P5 flow-cell end'), role('starter-derived PCR entry'),
                           role('tRNA-derived insert'), INDEX, role('P7 flow-cell end')),
                          name=TITLE)

def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=f'{CITATION}; doi:{DOI}; Methods, TGIRT and sequencing')
    trna = graph.input('trna', 'deacylated mature tRNA', MaterialKind.BIOLOGICAL,
                       state=molecule('deacylated tRNA', (role('structured tRNA'),
                                                         Segment('free 3-prime CCA', 'CCA', 'rna'))), primary=True)
    alkb = graph.input('alkb', 'wild-type AlkB plus AlkB-D135S', MaterialKind.CATALYST)
    demethylated = graph.product('demethylated_trna', 'demethylated tRNA', MaterialKind.MOLECULAR,
        state=molecule('demethylated tRNA', (role('structured tRNA with relieved RT blocks'),
                                           Segment('free 3-prime CCA', 'CCA', 'rna'))))
    graph.react('demethylate', 'Demethylate m1A, m3C and m1G', StepKind.CHEMICAL_CONVERSION,
                ((trna, InputRole.TEMPLATE), (alkb, InputRole.CATALYST)), (demethylated,))
    starter = graph.input('starter_duplex', 'blocked RNA / 32P-labeled DNA starter', MaterialKind.OLIGO,
                           state=starter_scene())
    tgirt = graph.input('tgirt', 'GsI-IIC MalE rigid fusion RT', MaterialKind.CATALYST)
    intermediate = graph.product('extended_starter', 'extended DNA with RNA starter still paired',
                                  MaterialKind.MOLECULAR, state=starter_scene(extended=True))
    graph.react('template_switch', 'INFERRED — DNA reconstruction extended by source-reported TGIRT switching',
                StepKind.TEMPLATE_SWITCHING,
                ((demethylated, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (tgirt, InputRole.CATALYST)), (intermediate,))
    alkali = graph.input('alkali', 'NaOH; neutralization with HCl', MaterialKind.CONDITION)
    cdna = graph.product('ss_cdna', 'purified single-stranded DNA', MaterialKind.MOLECULAR,
                         state=molecule('INFERRED — minimal starter-derived DNA', dna_parts()))
    rna_waste = graph.product('hydrolyzed_rna', 'hydrolyzed RNA', MaterialKind.WASTE,
                              state=molecule('RNA hydrolysis products', (role('RNA fragments'),)))
    graph.react('rna_hydrolysis', 'Hydrolyze RNA; recover gel-purified cDNA', StepKind.CHEMICAL_CONVERSION,
                ((intermediate, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)), (cdna, rna_waste))
    graph.discard(rna_waste, 'RNA is hydrolyzed before cDNA gel purification and circularization.')
    circligase = graph.input('circligase', 'CircLigase II (Epicentre)', MaterialKind.CATALYST)
    circle_construct = Construct(dna_parts(), name='minimal single-stranded cDNA')
    circle = graph.product('circular_cdna', 'circular single-stranded cDNA', MaterialKind.MOLECULAR,
                            state=closed_circle('INFERRED — minimal cDNA circle', dna_parts()))
    graph.react('circularize', 'INFERRED — minimal DNA assembly; source-reported CircLigase II circularization',
                StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known('single-strand intramolecular ligation',
                    statement='The 2015 Methods circularize purified single-stranded cDNA with CircLigase II.'))
    forward = graph.input('pcr_forward', 'source multiplex PCR primer', MaterialKind.OLIGO,
                           state=molecule('multiplex PCR primer', (DM_FORWARD,)))
    reverse = graph.input('pcr_reverse', 'source indexed PCR primer', MaterialKind.OLIGO,
                           state=molecule('indexed PCR primer; selected barcode unassigned', REVERSE_PARTS))
    polymerase = graph.input('pcr_polymerase', 'Phusion-HF (Thermo Scientific)', MaterialKind.CATALYST)
    library = graph.product('library', 'indexed DM-tRNA library', MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react('pcr', 'Amplify circular cDNA with source-published primers', StepKind.PCR,
                ((circle, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(circle_construct,
                    reason='Circular-to-linear PCR extension and selected barcode are not assembled in this model.'))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(tuple(graph.steps) + tuple(graph.materials), graph.authority)
    return graph.require_publishable()
