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

PAGE_FILE='user_cloning.html'
TITLE='USER cloning — illustrative uracil-excision example'
SOURCE='Bitinaite et al.2007, doi10.1093/nar/gkm041; published abstract, Introduction and Results USER enzyme composition/activity. Prepared vector entry only; its nicking/restriction preparation is outside this example.'
LEFT='ACGTAT'
RIGHT='AGCTAT'
A='CGATGACCTG'
B='TCAGGATCGA'

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    insert=CohesiveFragment('toy insert',LEFT,(base('insert body',A),),RIGHT,'3-prime','3-prime',inferred=True)
    vector=CohesiveFragment('prepared toy recipient',RIGHT,(base('recipient body',B),),LEFT,'3-prime','3-prime',inferred=True)
    canonical=LEFT+A+RIGHT
    top=LEFT[:-1]+'U'+A+RIGHT
    bot0=revcomp(canonical);bottom=bot0[:len(RIGHT)-1]+'U'+bot0[len(RIGHT):]
    source=g.input('insert_source','invented uracil-bearing duplex; U denotes dU',K.MOLECULAR,primary=True,state=duplex(Construct((base('toy uracil strand',top),)),bottom),inferred=True,source=EXAMPLE)
    recipient=g.input('recipient','invented already prepared recipient; not a replicating vector',K.MOLECULAR,primary=True,state=vector,chemical_state=vector,inferred=True,source=EXAMPLE)
    user=g.input('user','USER: uracil DNA glycosylase and endonuclease VIII',K.CATALYST)
    geometry=StrandRetention('insert_source','insert_ready',top,bottom,(len(LEFT),len(top)),(len(RIGHT),len(bottom)),insert,'Chosen uracil positions and retained strand intervals belong only to this invented example. Base removal plus backbone cleavage and terminal-fragment dissociation are separate chemical events; end-group microchemistry is outside this base-retention check.')
    ready=g.product('insert_ready','INFERRED — toy insert with 3-prime extensions',K.MOLECULAR,state=insert,chemical_state=insert)
    excised=g.product('excised','removed uracil-containing terminal intervals',K.MOLECULAR,state=Drawing((Row(chunks=[("Removed bases, abasic intermediates and short terminal fragments; exact chemical species unassigned.",None,False)]),)))
    g.react('excision','INFERRED — uracil excision, backbone cleavage and terminal-fragment release',S.DIGESTION,[(source,R.TEMPLATE),(user,R.CATALYST)],[ready,excised],geometry=geometry,relations=(M.STRAND_RETENTION,))
    g.mark_terminal(Terminal(excised.id,Disposition.UNRESOLVED,'Removed-base/abasic/short-fragment chemistry and physical cleanup are not expanded. These strands summarize removed reference intervals, not intact final chemical species.'))
    assembly=CircularCohesiveAssembly((('insert_ready',insert),('recipient',vector)),'Invented compatible 3-prime tails; actual recipient preparation remains outside scope.')
    con=assembly.product('toy reference circle')
    hybrid=g.product('hybrid','INFERRED — unsealed toy assembly',K.MOLECULAR,state=Drawing(tuple(assembly.rows())),chemical_state=con)
    g.react('anneal','INFERRED — anneal compatible extensions without in-vitro ligase',S.HYBRIDIZATION,[(ready,R.TEMPLATE),(recipient,R.STRUCTURAL)],[hybrid],geometry=assembly,relations=(M.COHESIVE_ASSEMBLY,))
    host=g.input('host','source transformation/recovery boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended host-repaired reference',K.MOLECULAR,state=closed_circle('synthetic repaired reference',tuple(con)),chemical_state=con)
    g.react('recover','INFERRED — host repair boundary',S.MAGIC,[(hybrid,R.TEMPLATE),(host,R.TEMPLATE)],[final],reason='Intracellular repair and clone selection are unexpanded. The short example has no replicon and is not a claim of recoverable clones.')
    g.constructed(final,'Illustrative circular-DNA reference; no NGS library or verified clone.')
    g.declare_reviewed_scope(('insert_source','recipient','user','excision','anneal','host','recover','final'),SOURCE+' Exact DNA and boundaries are invented.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_user_cloning.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/gkm041">Bitinaite et al.2007</a>. INFERRED — short invented DNA and exact uracil positions. U represents DNA deoxyuridine, not an RNA substrate. The recipient is already prepared; its actual source preparation is not replaced by PCR.'),
        WorkflowSection('Uracil-bearing insert → compatible extensions → repair boundary'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory toy relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
