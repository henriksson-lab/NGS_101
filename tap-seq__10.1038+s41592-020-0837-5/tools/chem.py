"""Bounded, source-local molecular workflow; not a complete protocol review."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition, SourceItem
from pcr import PCRGeometry
from page import head, info
PAGE_FILE='chem.html'

def unknown_template(label):
    return Construct([Segment(label,'X'*30,placeholder=True)])

def clean(g,key,material,label,selection,kind=K.MOLECULAR):
    sel=g.input(key+'_selection',selection,K.SELECTION)
    kept=g.product(key+'_kept',label,kind)
    waste=g.product(key+'_waste','Separated material from '+key,K.WASTE)
    g.react(key,label,S.CLEANUP,((material,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed only at this source-described physical separation; recovery is not predicted')
    return kept

def pcr(g,key,template,forward,reverse,label,conditions):
    f=g.input(key+'_forward',forward,K.OLIGO,state=primer_state(forward))
    r=g.input(key+'_reverse',reverse,K.OLIGO,state=primer_state(reverse))
    cat=g.input(key+'_enzyme',conditions,K.CATALYST)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,S.PCR,((template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(out,),
        geometry=PCRGeometry.unavailable(unknown_template('unassigned target sequence'),
            reason='This capped model does not assign the study-specific target substrate and complete primer-pool landing geometry.'))
    return out

def aliquot(g,key,material,label):
    selected=g.product(key+'_selected',label,K.MOLECULAR)
    remainder=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((material,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'The source specifies an aliquot; remaining reaction material is not claimed discarded.'))
    return selected
SOURCE='https://doi.org/10.1038/s41592-020-0837-5; Online Methods, Cell barcoding through Illumina sequencing'
TITLE='TAP-seq — 10x3primev2/CROP-seq branch'
LIMIT='Model begins at purified GEM-RT cDNA, preserving source cell identity and UMI. Historical bead sequence, target-panel landing geometry and exact final-strand bases remain unassigned. CROP-seq and gene targets share the targeted library.'

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cdna=g.input('cdna','35µl source10x3primev2 GEM-RT cDNA; RT template-switch oligo omitted upstream',K.MOLECULAR,primary=True)
    first=pcr(g,'outer',cdna,'Gene-specific outer pool plus8fold-excess CROPouter','PartialRead1 universal bead-handle primer','Outer targetedPCR products','KAPA HiFi;11cycles standard branch')
    first=clean(g,'outer_cleanup',first,'Purified outerPCR products','1.5×AMPureXP; two80%ethanol washes;30µl elution')
    second=pcr(g,'inner',aliquot(g,'outer_transfer',first,'10ng outerPCR product'),
        'Gene-specific inner pool plus8fold-excess CROPinner','PartialRead1 universal bead-handle primer','Semi-nested targetedPCR products','KAPA HiFi;8cycles')
    second=clean(g,'inner_cleanup',second,'Purified nestedPCR products','1.5×AMPureXP; two80%ethanol washes;30µl elution')
    third=pcr(g,'adapters',aliquot(g,'inner_transfer',second,'10ng nestedPCR product'),
        'Targeted10X primer','Selected Illumina N7XX reverse index primer','Illumina-adapted targetedPCR products','KAPA HiFi;8cycles')
    lib=clean(g,'adapter_cleanup',third,'TAP-seq gene/CROP-seq pooled targeted library','1.5×AMPureXP; two80%ethanol washes;30µl elution',K.LIBRARY)
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2','Index 1'),read_lengths={'Read 1':26,'Read 2':58,'Index 1':8},
        limitation='Main Methods: Read1 cell barcode+UMI26cycles; Read2 cDNA approximately58cycles; i7 sample index8cycles. Read2=58 is the representative approximate source schedule, not a fixed mandatory run length. Exact final strand geometry deferred.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def detail():return ''

REVIEWED=('cdna', 'outer_forward', 'outer_reverse', 'outer_enzyme', 'outer_product', 'outer_cleanup_selection', 'outer_cleanup_kept', 'outer_cleanup_waste', 'outer_transfer_selected', 'outer_transfer_remainder', 'inner_forward', 'inner_reverse', 'inner_enzyme', 'inner_product', 'inner_cleanup_selection', 'inner_cleanup_kept', 'inner_cleanup_waste', 'inner_transfer_selected', 'inner_transfer_remainder', 'adapters_forward', 'adapters_reverse', 'adapters_enzyme', 'adapters_product', 'adapter_cleanup_selection', 'adapter_cleanup_kept', 'adapter_cleanup_waste', 'outer', 'outer_cleanup', 'outer_transfer', 'inner', 'inner_cleanup', 'inner_transfer', 'adapters', 'adapter_cleanup')
def extra_checks(c,g):
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name
    original=source/SOURCE_FILE
    if original.exists():
        text=''.join(original.read_text().split())
        for oligo in TRANSCRIBED:
            c('Verified original transcription '+oligo,oligo in text)
    else:c.skip('Original transcription checks','External original unavailable')
    source_checks(c,g)


import illumina as il
import nextera as nx
PARTIAL_READ1='CTACACGACGCTCTTCCGATCT'
CROP_OUTER='AATGGACTATCATATGCTTACCGT'
CROP_INNER=nx.ADAPTOR_S7+'TGTGGAAAGGACGAAACACC'
TARGETED10X=il.P5+'TCTTTCCCTACACGACGCTCTTCCG'
SOURCE_FILE='s41592-020-0837-5_supp_EMS118361-supplement-Table_S6.xlsx.txt'
TRANSCRIBED=(PARTIAL_READ1,CROP_OUTER,CROP_INNER,TARGETED10X[:-3]+'C*C*G')
def primer_state(label):
    if label=='PartialRead1 universal bead-handle primer':return Construct([Segment('PartialRead1',PARTIAL_READ1)])
    if label=='Targeted10X primer':
        sc=Scene();sc.strand('primer',[Segment('Targeted10X',TARGETED10X)],label='Targeted10X; terminal C*C*G phosphorothioate linkages')
        return sc
    if label=='Selected Illumina N7XX reverse index primer':return Construct([
        Segment('P7',il.P7),Segment('sample index','N'*8,'idx',placeholder=True,feature=feature('tap.sample','sample_index','unknown')),
        Segment('N7 primer tail','GTCTCGTGGGCTCGG')])
    return None

def source_checks(c,g):
    c('Targeted guide transcript remains in same outerPCR pool','CROPouter' in g.materials['outer_forward'].label,True)
    c('Sample index role is structured',primer_state('Selected Illumina N7XX reverse index primer').get('sample index').feature.role,'sample_index')

def render_page():
    return '\n'.join([head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_tap_seq.html">Research notes</a></p>',
        info(LIMIT),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['adapter_cleanup_kept']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    c.report(exit_on_fail=True)
