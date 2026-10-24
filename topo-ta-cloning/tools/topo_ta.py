"""MAN0000047 RevisionA.0 pCR2.1-TOPO branch; invented neutral end model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Segment
from assembly_geometry import DuplexEndState,ActivatedEndTransfer
from cohesive_assembly import CohesiveFragment
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle,unsealed_circular_duplex
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='topo_ta.html'
TITLE='TOPO TA cloning — illustrative pCR2.1-TOPO end chemistry'
SOURCE='Thermo Fisher/Invitrogen TOPO TA Cloning Kit, MAN0000047 RevisionA.0, part250184; cover, Description of the system p9 and Produce PCR products p10. pCR2.1-TOPO subcloning branch only; no blunt-TOPO or Gateway-kit chemistry inherited.'
EXAMPLE='Invented short insert and recipient DNA bodies model the selected kit end chemistry. They are not the pCR2.1 sequence, a replicating vector, an experimental protocol instance or a validated construct.'

def body(name,seq):
    return Segment(name,seq,inferred=True)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    # The source-provided activated vector is the entry point: activation is not fabricated.
    vector=DuplexEndState(CohesiveFragment('toy recipient','A',(body('recipient body','TCAGGATCGATACG'),),'T','3-prime','3-prime',inferred=True),'unassigned','phosphotyrosyl','unassigned','phosphotyrosyl')
    insert=DuplexEndState(CohesiveFragment('toy A-tailed insert','T',(body('insert body','CGATGACCTGATAC'),),'A','3-prime','3-prime',inferred=True),'OH','OH','OH','OH')
    vm=g.input('vector','toy end model of source topoisomerase-activated recipient',K.MOLECULAR,primary=True,state=vector,chemical_state=vector,inferred=True,source=EXAMPLE)
    im=g.input('insert','toy supplied A-tailed PCR-product class; no 5-prime phosphates',K.MOLECULAR,primary=True,state=insert,chemical_state=insert,inferred=True,source=EXAMPLE)
    geom=ActivatedEndTransfer('vector','insert',vector,insert,'Source 3-prime phosphotyrosyl donors react with incoming physical 5-prime hydroxyl ends. This example checks declarations and reference order, not complete vector closure.')
    con=geom.reference('toy joined reference')
    joined=g.product('joined','INFERRED — intended transfer product; residual ends unassigned',K.MOLECULAR,state=unsealed_circular_duplex(tuple(con),statement='Activated vector donor ends accept insert 5-prime hydroxyls. Complete strand paths, remaining closure and exact nick positions are not assigned; no covalently closed reference strand is drawn.'),chemical_state=con)
    released=g.product('released','released topoisomerase class',K.MOLECULAR,state=Drawing((Row(chunks=[('Released enzyme; DNA-adduct/product distribution unassigned.',None,False)]),)))
    g.react('transfer','INFERRED — activated end transfer without added ligase',S.ADAPTER_ATTACHMENT,[(vm,R.TEMPLATE),(im,R.STRUCTURAL)],[joined,released],geometry=geom,relations=(M.ACTIVATED_END_TRANSFER,))
    g.mark_terminal(Terminal(released.id,Disposition.UNRESOLVED,'Released-enzyme handling and residual adducts are outside scope; no cleanup or complete conversion is claimed.'))
    host=g.input('host','source cell recovery and clone-analysis boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended repaired circular reference',K.MOLECULAR,state=closed_circle('intended toy reference',tuple(con)),chemical_state=con)
    g.react('recover','INFERRED — host recovery boundary',S.MAGIC,[(joined,R.TEMPLATE),(host,R.TEMPLATE)],[final],reason='Cellular repair and selection are not reconstructed. The toy lacks a replicon; no actual transformant, complete genotype or experimental success is claimed.')
    g.constructed(final,'Illustrative DNA reference only; source vector map, alternate insert orientation, PCR/terminal heterogeneity and NGS readiness remain outside this scope.')
    g.declare_reviewed_scope(('vector','insert','transfer','host','recover','final'),SOURCE+' Invented DNA bodies and chosen orientation are distinct from source end chemistry.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_topo_ta_cloning.html">Research notes</a></p>',
        info('<a href="https://assets.thermofisher.com/TFS-Assets/LSG/manuals/topota_man.pdf">MAN0000047 RevisionA.0</a>, printed pp9–10. INFERRED — invented DNA bodies and one illustrated orientation. The model starts with activated recipient and A-tailed product classes; it does not reconstruct vector activation or PCR.'),
        WorkflowSection('Activated 3-prime donors + 5-prime-OH insert → transfer → recovery boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory declared end-transfer relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/'topo-ta-cloning'/'topota_man.pdf.txt','Official MAN0000047 RevisionA.0')
    if have(c,source,label='selected vendor revision'):
        text=' '.join(source.path.read_text().split())
        c('source revision is the selected manual','Publication Number MAN0000047 Revision A.0' in text,True)
        c('source excludes phosphorylated PCR primers','Do not add 5´ phosphates to your primers for PCR' in text,True)
    c.report()
