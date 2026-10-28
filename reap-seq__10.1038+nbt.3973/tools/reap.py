"""Bounded REAP-seq protein-library model from the published supplement.

This starts after the source's 0.6x SPRI split. It is deliberately not a complete
cell-to-two-library graph. The inaccessible main Methods and the RNA branch remain
source-review obligations, rather than being filled from a modern 10x recipe.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, panel
from protocol_graph import (ProtocolGraph, InputRole as R, MaterialKind as K,
    StepKind as S, SequencingPlan, SourceItem)
from pcr import PCRGeometry
from page import head, info

PAGE_FILE = 'reap.html'
SOURCE = 'https://doi.org/10.1038/nbt.3973; published MOESM5, Discussion/Figs. 2-3/Tables 6-7'
S517 = 'GCGTAAGA'
AB_MOD5 = '/5AmMC6/'
AB_PRIMER = nx.READ1_PRIMER
R1_PRIMER = revcomp(il.TRUSEQ_READ2)
P5_INDEX_PRIMER = il.P5 + S517 + nx.READ1_PRIMER[:20]
P7_AS_WRITTEN = il.P7_RC  # Supplementary Table 6; do not silently relabel as canonical forward P7.
P5_AS_WRITTEN = il.P5

# Source Fig. 2: 25 dA, reduced to 24 when the preceding antibody barcode ends in A.
def antibody_tag(barcode):
    if len(barcode) != 8 or set(barcode)-set('ACGT'):
        raise ValueError('Antibody barcode must be eight concrete DNA bases')
    return AB_PRIMER + barcode + 'A'*(24 if barcode.endswith('A') else 25)

# Table 7 specifies order and read directions; resolving the Table 6 P7 PCR
# primer orientation against the inaccessible main Methods remains deferred.
# Accordingly all reconstructed strand coordinates default to inferred.
def protein_library():
    def s(name, seq, tag=None, **kw):
        return Segment(name, seq, tag, inferred=True, **kw)
    return Construct([
        s('P5', il.P5, 'p5'),
        s('S517 sample index', S517, 'i5', feature=feature('reap.sample','sample_index','fixed')),
        s('Nextera Read 1', AB_PRIMER, 'r1'),
        s('antibody barcode', 'N'*8, 'cbc', placeholder=True,
          feature=feature('reap.antibody','feature_barcode','unknown')),
        s('poly(dA)', 'A'*25),
        s('UMI', 'N'*10, 'umi', placeholder=True, feature=feature('reap.umi','umi','random')),
        s('TruSeq Read 2 reverse complement', R1_PRIMER, 'r2'),
        s('10x v1 cell barcode', 'N'*14, 'cbc', placeholder=True,
          feature=feature('reap.cell','cell_barcode','unknown')),
        s('P7 reverse complement', il.P7_RC, 'p7'),
    ], name='INFERRED — REAP protein strand coordinates from published Table 7')

PRIMERS = (
    sp.custom('Read 1','Custom TruSeq Read 2 reverse-complement primer', R1_PRIMER, SOURCE),
    sp.custom('Index 1 (i7)','Custom Nextera Read 1 primer', AB_PRIMER, SOURCE),
    sp.custom('Read 2','Standard TruSeq Read 2',il.TRUSEQ_READ2,SOURCE,
              note='Table 7 calls this standard Read 2; identity derived from the represented opposite-strand site.'),
)
I5 = sp.read_spec('Index 2 (i5)',8,'Sample index; standard primer and instrument orientation not resolved in this pass','demultiplex')
LENGTHS = {'Read 1':26,'Index 1 (i7)':8,'Read 2':10}

def protocol_graph():
    g = ProtocolGraph('REAP-seq — partial protein-library branch', authority=SOURCE)
    protein = g.input('protein_fraction','Short antibody-derived dsDNA in 0.6× SPRI supernatant (~155 bp), with free ssDNA',K.BIOLOGICAL,primary=True)
    exo = g.input('exonuclease','Exonuclease I',K.CATALYST)
    cleaned = g.product('exo_treated','dsDNA protein reporters after excess ssDNA digestion',K.MOLECULAR)
    g.react('exo_digest','Digest free ssDNA to reduce barcode crosstalk',S.DIGESTION,
            ((protein,R.TEMPLATE),(exo,R.CATALYST)),(cleaned,))
    p5 = g.input('p5_index','P5 / P5-Sample_Index-Part_Rd1 primers (Table 6)',K.OLIGO)
    p7 = g.input('p7','P7 primer as printed in Table 6; strand use unresolved',K.OLIGO)
    polymerase = g.input('polymerase','PCR polymerase; exact formulation/conditions outside inspected supplement scope',K.CATALYST)
    lib = g.product('protein_library','Indexed REAP protein library; strand reconstruction inferred',K.LIBRARY,
                    state=Scene.duplex(list(protein_library())),inferred=True,
                    source=SOURCE+'; Table 7 arrangement reconstructed; Table 6 P7 orientation still to resolve')
    g.react('protein_pcr','Amplify protein reporters with sample index and flow-cell adapters',S.PCR,
            ((cleaned,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(polymerase,R.CATALYST)),(lib,),
            geometry=PCRGeometry.unavailable(protein_library(), reason='Main Methods inaccessible; Table 6 prints a P7 reverse-complement sequence. Exact PCR primer polarity and amplification stages await reconciliation. No source-error claim is made.'))
    g.sequence(lib,SequencingPlan(library=protein_library(),primers=PRIMERS,
        acquired_roles=tuple(LENGTHS),read_lengths=LENGTHS,
        read_uses={'Read 1':'include','Index 1 (i7)':'include','Read 2':'include'},
        unplaced_reads=(I5,),
        limitation='This is only the protein branch. Coordinates are inferred from published Table 7; i5 standard priming orientation and upstream PCR primer polarity remain unresolved.'))
    g.declare_reviewed_scope((
        'protein_fraction','exonuclease','exo_treated','exo_digest','p5_index','p7',
        'polymerase','protein_library','protein_pcr'), SOURCE)
    g.add_source_item(SourceItem('cell_to_libraries','Cell labeling, droplet synthesis, size split and the complete RNA branch',
        SOURCE+'; Discussion and Figs. 1-3 identify these obligations, not represented by this partial graph'))
    return g

def render_page():
    lib = protein_library()
    return '\n'.join([head('REAP-seq — partial protein model'),'<div class="wrap"><h1>REAP-seq</h1>',
        '<p class="research-notes"><a href="01_reap_seq.html">Research notes</a></p>',
        info('INFERRED — Protein strand coordinates reconstructed from published Supplementary Table 7. This partial graph begins at the SPRI supernatant; the cell-to-library and RNA branches remain unmodelled. Exact PCR primer polarity and i5 instrument orientation require further source review.'),
        panel(Scene.duplex(list(lib)).rows(),cls='long',caption='INFERRED — protein-library strand reconstruction from Table 7; PCR polarity unresolved.'),
        sp.section(lib,PRIMERS,read_lengths=LENGTHS,unplaced_reads=(I5,)),
        '</div>'])

def selftest():
    from checks import Check
    c = Check()
    c('Source-specific barcode-ending-A manufacturing exception', len(antibody_tag('ACGTACGA')),65)
    c('Other antibody barcodes keep 25 dA',len(antibody_tag('ACGTACGT')),66)
    c('Table 6 P7 transcription retained in printed orientation',P7_AS_WRITTEN,'ATCTCGTATGCCGTCTTCTGCTTG')
    c('Table 6 index primer transcription',P5_INDEX_PRIMER,
      'AATGATACGGCGACCACCGAGATCTACACGCGTAAGATCGTCGGCAGCGTCAGATGT')
    c('Represented sequencing sites land',sp.verify(protein_library(),PRIMERS,required_roles=tuple(p.role for p in PRIMERS)),[])
    spans = [(x.read,x.feature.role,x.cycle_start,x.cycle_end)
             for x in sp.feature_spans(protein_library(),PRIMERS,LENGTHS)]
    c('Table 7 identifier-cycle assignments',spans,[
        ('Read 1','cell_barcode',1,14),('Index 1 (i7)','feature_barcode',1,8),
        ('Read 2','umi',1,10)])
    c('Partial graph cannot claim publication completeness',bool(protocol_graph().publication_problems()),True)
    from storage import sources_root
    source = sources_root()/Path(__file__).resolve().parents[1].name/'nbt.3973_41587_2017_BFnbt3973_MOESM5_ESM.pdf.txt'
    if source.exists():
        text=source.read_text()
        for seq in (AB_PRIMER,R1_PRIMER,P7_AS_WRITTEN,S517):
            c('Published supplement transcription: '+seq,seq in text)
    else:
        c.skip('Published supplement transcription','External supplement unavailable')
    c.report()
