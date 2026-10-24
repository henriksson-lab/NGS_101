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

from assembly_geometry import BridgedAssembly
PAGE_FILE='lcr.html'
TITLE='Ligase cycling reaction — illustrative bridge-oligo example'
SOURCE='de Kok et al.2014 doi10.1021/sb4001992; defining PubMed abstract PMID24932563. Abstract establishes bridging oligos, thermostable ligase and repeated denaturation/annealing/ligation; detailed original preparation remains deferred.'
A='ACGTATCGATGACCTG'
B='AGCTATTCAGGATCGA'
ARM=6

def example():
    a=Construct((base('body A',A),),'part A')
    b=Construct((base('body B',B),),'part B')
    bridges=(('bridge_ab',revcomp(A[-ARM:]+B[:ARM])),('bridge_ba',revcomp(B[-ARM:]+A[:ARM])))
    return BridgedAssembly((('a',a),('b',b)),bridges,ARM,True,True,'Invented blunt parts and bridge arms; both DNA strands are declared 5-prime phosphorylated and 3-prime hydroxyl terminated for this illustration. This does not establish the original paper preparation or experimental annealing.')

def bridge_population(geom):
    sc=Scene()
    for identifier,sequence in geom.bridges:
        sc.strand(identifier,(base(identifier,sequence),),label='bridge association unassigned')
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    geom=example();fragments=[];bridges=[]
    for identifier,con in geom.fragments:
        fragments.append(g.input(identifier,'invented 5-prime phosphorylated blunt duplex '+identifier,K.MOLECULAR,primary=True,state=duplex(con),chemical_state=con,inferred=True,source=EXAMPLE))
    for identifier,sequence in geom.bridges:
        bridges.append(g.input(identifier,'invented nonincorporated bridging oligo',K.OLIGO,state=Construct((base(identifier,sequence),)),inferred=True,source=EXAMPLE))
    ligase=g.input('ligase','source thermostable DNA ligase',K.CATALYST)
    product=geom.product('toy circular reference')
    final=g.product('final','INFERRED — ligated synthetic circular reference',K.MOLECULAR,state=closed_circle('toy circular reference',tuple(product)),chemical_state=product)
    bridge_pool=g.product('remaining_bridges','bridging oligos; association/partition after cycling unassigned',K.MOLECULAR,state=bridge_population(geom))
    g.react('cycle','INFERRED — bridge-directed annealing and ligation within thermal cycling',S.LIGATION,[(fragments[0],R.TEMPLATE),(fragments[1],R.STRUCTURAL)]+[(m,R.STRUCTURAL) for m in bridges]+[(ligase,R.CATALYST)],[final,bridge_pool],geometry=geom,relations=(M.BRIDGED_ASSEMBLY,))
    g.mark_terminal(Terminal(bridge_pool.id,Disposition.UNRESOLVED,'Bridge DNA is excluded from the product reference. Its binding and removal after cycling are not assigned; no cleanup is inferred.'))
    g.constructed(final,'Only the declared DNA parts enter the circular reference; no bridge bases are covalently incorporated.')
    g.declare_reviewed_scope(('a','b','bridge_ab','bridge_ba','ligase','cycle','final'),SOURCE+' Exact DNA, arm lengths and prepared end chemistry are invented example preconditions.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_lcr.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1021/sb4001992">de Kok et al.2014</a>. INFERRED — invented blunt fragments, bridge arms and prepared end chemistry. Short arm lengths illustrate geometry, not a workable design.'),
        WorkflowSection('Bridge-supported seams and fragment-only circular product'),
        info('The graph summarizes the molecular relation within source thermal cycling; it does not claim that one cycle is sufficient or reconstruct the full original workflow.'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();claims=check_graph(protocol_graph())
    c('no contradictory toy molecular relation',[x for x in claims if x.status is Status.CONTRADICTORY],[])
    c('bridge relation checked independently',sum(x.question==M.BRIDGED_ASSEMBLY.value for x in claims),1)
    c.report()
