from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Row, Scene, Segment, revcomp
from cohesive_assembly import CohesiveFragment, CircularCohesiveAssembly
from assembly_geometry import StrandRetention
from protocol_graph import ProtocolGraph, Disposition, Terminal, InputRole as R, MaterialKind as K, StepKind as S, MolecularRelation as M
from protocol_model import Drawing, closed_circle
from protocol_views import WorkflowSection, assemble_page
from page import head,info

REPRESENTATION='illustrative_example'
EXAMPLE='Invented short, nonfunctional DNA parts demonstrate the declared molecular relationships. They were not used in the cited experiments and are not a working vector or validated laboratory design.'

def base(name,sequence):
    return Segment(name,sequence,inferred=True)

def duplex(top,bottom=None):
    topseq=top.top(); bottomseq=revcomp(topseq) if bottom is None else bottom
    sc=Scene();sc.strand('top',(base('top reference',topseq),))
    sc.anneal('bottom',(base('bottom reference',bottomseq),),to='top',pair=('bottom reference','top reference'))
    return sc

def remnants(geometry):
    count=sum(map(len,geometry.removed()))
    return Drawing((Row(chunks=[(f"Removed parent reference intervals: {count} strand positions. Cleavage/degradation products are unexpanded; no intact excised strand is implied.",None,False)]),))

PAGE_FILE='gibson_assembly.html'
TITLE='Gibson assembly — illustrative DNA example'
SOURCE='Gibson et al.2009; defining abstract and published Supplementary Text: separate 5-prime isothermal method from 3-prime alternatives.'
DOI='10.1038/nmeth.1318'

# Chosen only for a small symbolic example, not an experimental design.
LEFT='ACGTAT'
RIGHT='AGCTAT'
A='CGATGACCTG'
B='TCAGGATCGA'

def example_fragments():
    return (CohesiveFragment('part A',LEFT,(base('body A',A),),RIGHT,'3-prime','3-prime',inferred=True),
            CohesiveFragment('part B',RIGHT,(base('body B',B),),LEFT,'3-prime','3-prime',inferred=True))

def precursor(fragment):
    core=''.join(s.top for s in fragment.core);top=fragment.left+core+fragment.right
    parent=Construct((base('invented overlap fragment',top),),'overlapping toy duplex')
    left=len(fragment.left);right=len(fragment.right)
    if fragment.left_polarity=='3-prime':
        return parent,(left,len(top)),(right,len(top))
    return parent,(0,len(top)-right),(0,len(top)-left)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    enzyme=g.input('enzyme','source 5-prime exonuclease',K.CATALYST)
    prepared=[]
    for identifier,fragment in zip(('a','b'),example_fragments()):
        parent,top_span,bottom_span=precursor(fragment)
        source=g.input(identifier+'_source','invented duplex '+identifier,K.MOLECULAR,primary=True,state=duplex(parent),inferred=True,source=EXAMPLE)
        geom=StrandRetention(source.id,identifier+'_ready',parent.top(),revcomp(parent.top()),top_span,bottom_span,fragment,'Explicit synthetic strand intervals; not a prediction of enzyme stopping distributions.')
        ready=g.product(identifier+'_ready','INFERRED — processed toy part '+identifier,K.MOLECULAR,state=fragment,chemical_state=fragment)
        removed=g.product(identifier+'_removed','toy removed strand intervals',K.MOLECULAR,state=remnants(geom))
        g.react(identifier+'_process','INFERRED — '+'5-prime resection; exact stopping positions are invented',S.DIGESTION,[(source,R.TEMPLATE),(enzyme,R.CATALYST)],[ready,removed],geometry=geom,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Removed parent bases remain represented by an unresolved product class. Molecular cleavage/degradation products and physical removal are not expanded.'))
        prepared.append(ready)
    geometry=CircularCohesiveAssembly(tuple((m.id,f) for m,f in zip(prepared,example_fragments())),'Declared toy order and complementary end polarities; only reference sequence and seam pairing are checked.')
    product=geometry.product('toy circular reference')
    joint=g.product('joint','INFERRED — paired toy intermediate; two nicks per seam',K.MOLECULAR,state=Drawing(tuple(geometry.rows())),chemical_state=product)
    g.react('anneal','INFERRED — anneal the declared compatible toy ends',S.HYBRIDIZATION,[(prepared[0],R.TEMPLATE),(prepared[1],R.STRUCTURAL)],[joint],geometry=geometry,relations=(M.COHESIVE_ASSEMBLY,))

    polymerase=g.input('polymerase','source DNA polymerase',K.CATALYST)
    repaired=g.product('repaired','INFERRED — polymerase-repair boundary; nicks remain',K.MOLECULAR,state=Drawing(tuple(geometry.rows())),chemical_state=product)
    g.react('repair','INFERRED — polymerase repair; zero gaps in this exact-end toy',S.EXTENSION,[(joint,R.TEMPLATE),(polymerase,R.CATALYST)],[repaired],relations=(M.STRAND_ORDER,),reason='The source mixture includes polymerase. This limiting toy retains exactly matching overhangs, so no gap or copied nucleotide is invented. Variable resection and gapped assemblies need separate explicit geometry.')
    ligase=g.input('ligase','source DNA ligase',K.CATALYST)
    final=g.product('final','INFERRED — sealed synthetic circular DNA',K.MOLECULAR,state=closed_circle('toy circular DNA',tuple(product)),chemical_state=product)
    g.react('seal','INFERRED — seal remaining nicks',S.NICK_SEALING,[(repaired,R.TEMPLATE),(ligase,R.CATALYST)],[final],relations=(M.STRAND_ORDER,))

    g.constructed(final,'Illustrative constructed-DNA reference only, not an experimental clone or NGS library.')
    g.declare_reviewed_scope(('enzyme','a_source','a_process','b_source','b_process','anneal','polymerase','repair','ligase','seal','final'),SOURCE+' Example base values and exact processing boundaries are invented and separate from the source method.')
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_gibson_assembly.html">Research notes</a></p>',
        info('<a href="https://doi.org/'+DOI+'">Defining method</a>. INFERRED — invented bases and exact processing boundaries; known method chemistry is described separately in the notes.'),
        WorkflowSection('Illustrative molecular steps'),
        info('The panels separate molecular relations, not separate incubations. Source alternatives, experimental conditions, yields and complete vector construction are outside this example.'),'</div>'))


def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();claims=check_graph(protocol_graph())
    c('no contradictory toy molecular claims',[x for x in claims if x.status is Status.CONTRADICTORY],[])
    c('retained strand and end-assembly relations are independently checked',sum(x.question==M.STRAND_RETENTION.value for x in claims),2)
    c.report()
