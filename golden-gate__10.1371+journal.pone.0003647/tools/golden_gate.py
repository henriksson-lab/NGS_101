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

PAGE_FILE='golden_gate.html'
TITLE='Golden Gate — illustrative DNA example'
SOURCE='Engler et al.2008; published abstract and Results basic cloning strategy, Figure1 caption.'
DOI='10.1371/journal.pone.0003647'

# Chosen only for a small symbolic example, not an experimental design.
LEFT='ACGTAT'
RIGHT='AGCTAT'
A='CGATGACCTG'
B='TCAGGATCGA'

def example_fragments():
    return (CohesiveFragment('part A',LEFT,(base('body A',A),),RIGHT,'5-prime','5-prime',inferred=True),
            CohesiveFragment('part B',RIGHT,(base('body B',B),),LEFT,'5-prime','5-prime',inferred=True))

from plasmid import BsaI,Plasmid,digest
# BsaI makes four-base ends; these invented fusion choices replace the six-base overlap toy above.
LEFT='AGTC'
RIGHT='CTGA'

def precursor(fragment):
    top=BsaI.site+'A'+fragment.left+''.join(s.top for s in fragment.core)+fragment.right+'T'+revcomp(BsaI.site)
    parent=Construct((base('invented BsaI-flanked part',top),),'BsaI-flanked toy')
    cuts=digest(Plasmid('toy',top,circular=False),BsaI)
    if len(cuts)!=2 or cuts[0].strand!=1 or cuts[1].strand!=-1:
        raise ValueError('toy requires the declared inward-cutting BsaI pair only')
    lo,hi=cuts
    return parent,(lo.top_cut,hi.top_cut),(len(top)-hi.bottom_cut,len(top)-lo.bottom_cut)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    enzyme=g.input('enzyme','BsaI',K.CATALYST)
    prepared=[]
    for identifier,fragment in zip(('a','b'),example_fragments()):
        parent,top_span,bottom_span=precursor(fragment)
        source=g.input(identifier+'_source','invented duplex '+identifier,K.MOLECULAR,primary=True,state=duplex(parent),inferred=True,source=EXAMPLE)
        geom=StrandRetention(source.id,identifier+'_ready',parent.top(),revcomp(parent.top()),top_span,bottom_span,fragment,'Explicit synthetic strand intervals; not a prediction of enzyme stopping distributions.')
        ready=g.product(identifier+'_ready','INFERRED — processed toy part '+identifier,K.MOLECULAR,state=fragment,chemical_state=fragment)
        removed=g.product(identifier+'_removed','toy removed strand intervals',K.MOLECULAR,state=remnants(geom))
        g.react(identifier+'_process','INFERRED — '+'BsaI cleavage within one-pot restriction–ligation',S.DIGESTION,[(source,R.TEMPLATE),(enzyme,R.CATALYST)],[ready,removed],geometry=geom,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Removed parent bases remain represented by an unresolved product class. Molecular cleavage/degradation products and physical removal are not expanded.'))
        prepared.append(ready)
    geometry=CircularCohesiveAssembly(tuple((m.id,f) for m,f in zip(prepared,example_fragments())),'Declared toy order and complementary end polarities; only reference sequence and seam pairing are checked.')
    product=geometry.product('toy circular reference')
    joint=g.product('joint','INFERRED — paired toy intermediate; two nicks per seam',K.MOLECULAR,state=Drawing(tuple(geometry.rows())),chemical_state=product)
    g.react('anneal','INFERRED — anneal the declared compatible toy ends',S.HYBRIDIZATION,[(prepared[0],R.TEMPLATE),(prepared[1],R.STRUCTURAL)],[joint],geometry=geometry,relations=(M.COHESIVE_ASSEMBLY,))

    ligase=g.input('ligase','source T4 DNA ligase',K.CATALYST)
    final=g.product('final','INFERRED — sealed synthetic circular DNA',K.MOLECULAR,state=closed_circle('toy circular DNA',tuple(product)),chemical_state=product)
    g.react('seal','INFERRED — ligate the toy junctions within the same mixture',S.NICK_SEALING,[(joint,R.TEMPLATE),(ligase,R.CATALYST)],[final],relations=(M.STRAND_ORDER,))
    if BsaI.site in product.top()+product.top()[:len(BsaI.site)-1] or revcomp(BsaI.site) in product.top()+product.top()[:len(BsaI.site)-1]:
        raise ValueError('intended toy product must lose both orientations of the recognition site')

    g.constructed(final,'Illustrative constructed-DNA reference only, not an experimental clone or NGS library.')
    g.declare_reviewed_scope(('enzyme','a_source','a_process','b_source','b_process','anneal','ligase','seal','final'),SOURCE+' Example base values and exact processing boundaries are invented and separate from the source method.')
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_golden_gate.html">Research notes</a></p>',
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
