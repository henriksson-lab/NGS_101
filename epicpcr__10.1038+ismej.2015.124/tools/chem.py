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
    f=g.input(key+'_forward',forward,K.OLIGO)
    r=g.input(key+'_reverse',reverse,K.OLIGO)
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
SOURCE='https://doi.org/10.1038/ismej.2015.124; Methods and Supplementary Methods/TablesS2,S4,S5'
TITLE='epicPCR — dsrB–16S branch'
LIMIT='Fusion linkage is source-declared, not a prediction of specificity or occupancy. Unfused products are blocked in reaction and removed only at physical cleanups. Full unknown dsrB/16S templates and sequencing-primer coordinates remain unassigned.'
F1='GTGTAGCAGTTACCGCA'
BRIDGE='GWATTACCGCGGCKGCTGTGCCTSAAYATGTGYGGYG'
NESTED='ACACGACGCTCTTCCGATCTYRYRVAGVATSGCGATRTCGGA'
BLOCK_R='TTTTTTTTTTGWATTACCGCGGCKGCTG'
BLOCK_F='TTTTTTTTTTCAGCMGCCGCGGTAATWC'

def blocker_scene():
    sc=Scene()
    sc.strand('R',[Segment('U519R-block10',BLOCK_R,placeholder=True)],mod3='C3 spacer',label='U519R-block10 (3′ extension blocked)')
    sc.strand('F',[Segment('U519F-block10',BLOCK_F,placeholder=True)],mod3='C3 spacer',label='U519F-block10 (3′ extension blocked)')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    beads=g.input('beads','Lysed environmental cells in source-prepared polyacrylamide hydrogel beads',K.BIOLOGICAL,primary=True)
    f=g.input('fusion_forward','dsrB-F1',K.OLIGO,state=Construct([Segment('dsrB-F1',F1)]))
    b=g.input('bridge','Limiting dsrB-R1_519R fusion primer:10nM',K.OLIGO,state=Construct([Segment('dsrB-R1_519R',BRIDGE,placeholder=True)]))
    r=g.input('fusion_reverse','Universal16S R2 primer from source selected panel',K.OLIGO)
    enzyme=g.input('fusion_enzyme','Source fusionPCR mix in ABILem90 emulsion; F1/R2 each1µM',K.CATALYST)
    fused=g.product('fused','Fusion reaction:dsrB–16S products plus unfused/off-target material',K.MOLECULAR)
    g.react('fusion','Compartmentalized three-primer overlap-fusionPCR',S.PCR,
        ((beads,R.TEMPLATE),(f,R.STRUCTURAL),(b,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(fused,),
        geometry=PCRGeometry.unavailable(unknown_template('dsrB and16S templates'),reason='Bridge sequence known, but study-specific genomic target sequences and full overlap-extension geometry are unassigned.'))
    recovered=clean(g,'emulsion_recovery',fused,'Recovered fusion-PCR DNA','Break emulsion, transfer aqueous fraction, AMPure cleanup')
    nested_f=g.input('nested_forward','i_dsrB-F3; includes YRYR sequence-diversity bases',K.OLIGO,state=Construct([Segment('i_dsrB-F3',NESTED,placeholder=True)]))
    nested_r=g.input('nested_reverse','i_E786R nested16S reverse primer',K.OLIGO)
    blocks=g.input('blockers','U519R-block10 and U519F-block10 with3prime C3 spacers;3.2µM each',K.OLIGO,state=blocker_scene())
    cat=g.input('nested_enzyme','Source nestedPCR mix; cycle count chosen using qPCR',K.CATALYST)
    nested=g.product('nested_product','Nested fusion-PCR reaction; unfused material not assumed absent',K.MOLECULAR)
    g.react('nested','NestedPCR with nonextendable blocking oligos',S.PCR,
        ((recovered,R.TEMPLATE),(nested_f,R.STRUCTURAL),(nested_r,R.STRUCTURAL),(blocks,R.STRUCTURAL),(cat,R.CATALYST)),(nested,),
        geometry=PCRGeometry.unavailable(unknown_template('fusion amplicon'),reason='C3 blocking chemistry is source-declared; this graph does not simulate competition, A addition, or selective amplification.'))
    purified=clean(g,'nested_cleanup',nested,'Cleaned nested fusion amplicons','Source postnested cleanup')
    final=pcr(g,'adapter_pcr',purified,'PE-PCR-F forward primer','Selected PE-PCR-XXX sample-index reverse primer','Indexed dsrB–16S amplicons','Source adapter-addition PCR')
    library=clean(g,'library_cleanup',final,'epicPCR dsrB–16S library','Source final library purification',K.LIBRARY)
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':250,'Read 2':250},
        limitation='Source MiSeq2×250; selected sample index and exact read-primer coordinates remain unassigned. YRYR is diversity sequence, not a UMI.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def detail():return panel(blocker_scene().rows(),caption='Source blocking oligos carry a nonextendable3′ C3 spacer')

REVIEWED=('beads', 'fusion_forward', 'bridge', 'fusion_reverse', 'fusion_enzyme', 'fused', 'emulsion_recovery_selection', 'emulsion_recovery_kept', 'emulsion_recovery_waste', 'nested_forward', 'nested_reverse', 'blockers', 'nested_enzyme', 'nested_product', 'nested_cleanup_selection', 'nested_cleanup_kept', 'nested_cleanup_waste', 'adapter_pcr_forward', 'adapter_pcr_reverse', 'adapter_pcr_enzyme', 'adapter_pcr_product', 'library_cleanup_selection', 'library_cleanup_kept', 'library_cleanup_waste', 'fusion', 'emulsion_recovery', 'nested', 'nested_cleanup', 'adapter_pcr', 'library_cleanup')
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


SOURCE_FILE='ismej.2015.124_supp_ismej2015124x1.pdf.txt'
TRANSCRIBED=(F1,BRIDGE,NESTED,BLOCK_R+'/3SpC3/',BLOCK_F+'/3SpC3/')
def source_checks(c,g):
    c('Blockers are present in nested reaction, not discarded beforehand','blockers' in tuple(x.material for x in g.steps['nested'].inputs),True)
    c('Blocking oligos have nonextendable ends',tuple(x.mod3 for x in blocker_scene().strands.values()),('C3 spacer','C3 spacer'))

def render_page():
    return '\n'.join([head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_epicpcr.html">Research notes</a></p>',
        info(LIMIT),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['library_cleanup_kept']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    c.report(exit_on_fail=True)
