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

PAGE_FILE='slice.html'
TITLE='SLiCE — illustrative extract-mediated assembly'
SOURCE='Zhang, Werling and Edelmann2012 doi10.1093/nar/gkr1288; published abstract, Introduction and Discussion generic end-homology/extract mechanism. Exact extract activities and molecular intermediates are not assigned.'
LEFT='ACGTAT'
RIGHT='AGCTAT'
A='CGATGACCTG'
B='TCAGGATCGA'

def references():
    a=Construct((base('left overlap',LEFT),base('body A',A),base('right overlap',RIGHT)),'toy part A')
    b=Construct((base('left overlap',RIGHT),base('body B',B),base('right overlap',LEFT)),'toy part B')
    # Explicit example overlap identity; this does not predict an extract reaction.
    if a.top()[:len(LEFT)]!=b.top()[-len(LEFT):] or a.top()[-len(RIGHT):]!=b.top()[:len(RIGHT)]:
        raise ValueError('invented overlap identities disagree')
    product=Construct((base('overlap 1',LEFT),base('body A',A),base('overlap 2',RIGHT),base('body B',B)),'intended toy reference')
    return a,b,product

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    a,b,con=references()
    ma=g.input('a','invented overlapping duplex A',K.MOLECULAR,primary=True,state=duplex(a),inferred=True,source=EXAMPLE)
    mb=g.input('b','invented overlapping duplex B',K.MOLECULAR,primary=True,state=duplex(b),inferred=True,source=EXAMPLE)
    extract=g.input('extract','source bacterial cell extract; activities unresolved',K.CATALYST)
    joined=g.product('joined','INFERRED — intended joined reference; chemical state unresolved',K.MOLECULAR,state=Drawing((Row(chunks=[('Extract product molecular topology, strand ends and covalent closure are unassigned; the intended reference order is conditional.',None,False)]),)))
    g.react('extract_reaction','INFERRED — extract-mediated assembly boundary',S.MAGIC,[(ma,R.TEMPLATE),(mb,R.TEMPLATE),(extract,R.CATALYST)],[joined],reason='Source establishes extract-mediated end-homology assembly. No nuclease direction, exact resection positions, added purified ligase or defined gap-repair recipe is assigned to the extract. The small example only declares an intended overlap-preserving reference order.')
    host=g.input('host','source transformation and recovery boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended recovered circular reference',K.MOLECULAR,state=closed_circle('intended synthetic reference',tuple(con)),chemical_state=con)
    g.react('recover','INFERRED — host recovery boundary',S.MAGIC,[(joined,R.TEMPLATE),(host,R.TEMPLATE)],[final],reason='Neither intracellular repair nor extraction-versus-host contributions are resolved. The toy has no replicon, so no actual transformant or experimental success is claimed.')
    g.constructed(final,'Conditional intended example sequence only; no extract activity, experimental recovery or NGS readiness is established.')
    g.declare_reviewed_scope(('a','b','extract','extract_reaction','host','recover','final'),SOURCE+' Exact short DNA and intended example order are invented.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_slice.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/gkr1288">Zhang et al.2012</a>. INFERRED — invented short parts and intended reference order. The extract and host boundaries retain their uncertainty; no exact chemical intermediate is asserted.'),
        WorkflowSection('Overlapping parts → extract boundary → intended reference'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory declared molecular relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
