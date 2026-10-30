from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from protocol_graph import Terminal,Disposition

def save(g,key,x,label):
    take=g.product(key+'_aliquot',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Remaining preparation; later use unassigned',K.MOLECULAR)
    g.react(key,'Physical aliquot: '+label,S.GENERIC,((x,R.TEMPLATE),),(take,rest))
    g.terminals[rest.id]=Terminal(rest.id,Disposition.STORED,'Remainder retained explicitly; later fate unassigned.')
    return take

def render_page():
    g=protocol_graph()
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(g,ENDPOINT))+'</main>'

TITLE='DNA-barcoded pMHC multimers — published reporter components'
SOURCE='Published MOESM2 FigS1/Table1;MOESM3 Table2;MOESM4 Table3'
LIMIT='Private source-limited reporter graph. Defining full main unavailable in standard retrieval. Published SI supports reporter assembly and sorted-cell barcode recovery at overview level;core extension/staining/recovery/PCR/cleanup and run order remain unreviewed. Do not borrow TetTCR or current vendor chemistry.'
ENDPOINT='library'
SOURCE_OLIGO_EXCEPTIONS={'oligo_a':'Published nbt.3662 MOESM2 Table1 originalPDF p11; extracted lines119–126 interpose row label A1 between GGTCAGCA and TCATTTCC. Whole BiotinA1 strand including6N/25nt barcode/16nt arm independently checked against rendered original by two readers2026-10-29. Narrow table-extraction exception;not publication error or absence audit.'}
A_FIXED='GAAGTTCCAGCCAGCGTCACAGTTT';A_BAR='CGAGGGCAATGGTTAACTGACACGT';A_TAIL='GGTCAGCATCATTTCC'
B_FIXED='CTGTGACTATGTGAGGCTTTCTCGA';B_BAR='GCCTGTAGTCCCACGCGATCTAAC A'.replace(' ','');B_TAIL='GGAAATGATGCTGACC'
def oligo_a():
    sc=Scene();sc.strand('a',[Segment('source primer region',A_FIXED),Segment('source UMI','N'*6,feature=feature('pmhc.a.umi','umi','random')),Segment('source25nt barcode',A_BAR,feature=feature('pmhc.a1','feature_barcode','fixed',note='Source reporter identity,not cell identity.')),Segment('annealing arm',A_TAIL)],mod5='Biotin',label='Published biotin-Oligo A1');return sc
def oligo_b():
    return Construct([Segment('source primer region',B_FIXED),Segment('source UMI','N'*6,feature=feature('pmhc.b.umi','umi','random')),Segment('source25nt barcode',B_BAR,feature=feature('pmhc.b1','feature_barcode','fixed')),Segment('annealing arm',B_TAIL)])
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    a=g.input('oligo_a','Published biotin OligoA1',K.OLIGO,state=oligo_a());b=g.input('oligo_b','Published OligoB1',K.OLIGO,state=oligo_b());back=g.input('backbone','Source multimerization backbone;exact assembly protocol unreviewed',K.MOLECULAR)
    reporter=g.product('reporter','DNA-barcoded multimerization backbone preparation;published caption says two identical barcodes per backbone,actual occupancy not predicted',K.MOLECULAR)
    g.react('reporter_assembly','Source DNA barcode attachment before pMHC co-attachment;full extension/copy/end geometry unreviewed',S.MAGIC,((a,R.TEMPLATE),(b,R.STRUCTURAL),(back,R.STRUCTURAL)),(reporter,),reason='Published SI FigS1 states assembly order;Table1/2 define components,not full molecular operations.')
    pmhc=g.input('pmhc','Source peptide-MHC molecules;UV exchange/manufacture unreviewed',K.MOLECULAR)
    mult=g.product('multimer','Source DNA-barcoded pMHC multimer preparation,unbound classes not declared removed',K.MOLECULAR)
    g.react('pmhc_assembly','Co-attach pMHC to source DNA-labelled backbone',S.GENERIC,((reporter,R.TEMPLATE),(pmhc,R.STRUCTURAL)),(mult,))
    cells=g.input('cells','Source T-cell population',K.BIOLOGICAL,primary=True)
    stain=g.product('stained','Source multimer-stained cells and unbound components;core workup unreviewed',K.MOLECULAR)
    g.react('stain','Published SI reported multimer staining;exact core order unreviewed',S.MAGIC,((cells,R.TEMPLATE),(mult,R.STRUCTURAL)),(stain,),reason='Only published SI experiment overview inspected;no ideal antigen recognition claim.')
    selected=g.product('selected','Reported sorted population with associated barcode DNA;core recovery chemistry unreviewed',K.MOLECULAR);other=g.product('other_cells','Other sorted/rejected fractions;experiment-specific later fate unassigned',K.MOLECULAR)
    g.react('sorting','Published SI FACS sorted-cell population selection',S.CLEANUP,((stain,R.TEMPLATE),),(selected,other));g.terminals[other.id]=Terminal(other.id,Disposition.UNRESOLVED,'Unselected populations not assigned library membership.')
    f=g.input('forward','Published F1 sample-ID/A-key primer;feature boundaries unreviewed',K.OLIGO,state=Construct([Segment('source F1','CCATCTCATCCCTGCGTGTCTCCGACTCAGCTGGGGGAAGTTCCAGCCAGCGTC')]))
    r=g.input('reverse','Published R1 P1-key primer',K.OLIGO,state=Construct([Segment('source R1','CCTCTCTATGGGCAGTCGGTGATCTGTGACTATGTGAGGCTTTC')]))
    lib=g.product('library','Reported barcode-library continuation;DNA copies,source-bound tags and handling unassembled',K.MOLECULAR)
    g.react('library_continuation','Unreviewed defining-source recovery/library operations',S.MAGIC,((selected,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),(lib,),reason='Table3 identifies primers,not reaction/cleanup order;main Methods unavailable. No inferred one-step PCR or invented thermal recipe.')
    end(g,lib,'Published Table3 identifies IonTorrent A-Key/P1-Key adapters;instrument/cycles/run-specific acquisition not inspected.')
    return incomplete(g)
def selftest():
    from checks import Check
    from storage import sources_root
    from chemdraw import revcomp
    c=Check();g=protocol_graph()
    c('Source annealing arms complementary',revcomp(A_TAIL),B_TAIL)
    c('Biotin retained on actual molecular strand',oligo_a().strands['a'].mod5,'Biotin')
    c('UMI distinguished from analyte barcode',oligo_a().strands['a'].segs[1].feature.role,'umi')
    c('No barcode interpreted as cell identity',oligo_b().segments[2].feature.role,'feature_barcode')
    c('Unselected cell fraction not silently library',g.terminals['other_cells'].disposition,Disposition.UNRESOLVED)
    c('Main unreviewed prevents publication',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nbt.3662_41587_2016_BFnbt3662_MOESM3_ESM.xlsx.txt'
    if p.exists():c('Whole B1 source transcription',oligo_b().top() in p.read_text(),True)
    else:c.skip('Original Table2','External original absent')
    c.report()
