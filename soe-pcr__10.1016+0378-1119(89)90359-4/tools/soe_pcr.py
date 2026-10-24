"""SOE-PCR: invented linear fusion with source-qualified stages and exact geometry."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,revcomp
from pcr import PCRGeometry,PCRPrimer
from assembly_geometry import pcr_product,OverlapExtension
from protocol_graph import ProtocolGraph,Disposition,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='soe_pcr.html'
TITLE='SOE-PCR — illustrative linear overlap fusion'
SOURCE='Horton et al.1989 doi10.1016/0378-1119(89)90359-4; defining PubMed abstract PMID2744488. Separate PCR fragments, complementary ends and mutual extension; detailed original methods deferred.'
REPRESENTATION='illustrative_example'
EXAMPLE='Invented short neutral DNA and primers, not used in the cited experiments. This is a reference-copying example, not a validated laboratory design or functional construct.'
LEFT='ACGTCAGTACGATCGTGACTACGA'
OVERLAP='GATCTACGAGTC'
RIGHT='TCGAGTACCTGACATGTCAGATCA'

def dna(name,sequence):
    return Construct((Segment(name,sequence,inferred=True),),name)

def duplex(con):
    return Scene.duplex(tuple(con),label=con.name)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation=REPRESENTATION,example_description=EXAMPLE)
    a=g.input('a','invented starting template A',K.MOLECULAR,primary=True,state=duplex(dna('body A',LEFT)),chemical_state=dna('body A',LEFT),inferred=True,source=EXAMPLE)
    b=g.input('b','invented starting template B',K.MOLECULAR,primary=True,state=duplex(dna('body B',RIGHT)),chemical_state=dna('body B',RIGHT),inferred=True,source=EXAMPLE)
    pol=g.input('polymerase','source DNA polymerase',K.CATALYST)
    primers={}
    for identifier,sequence in (('outer_f',LEFT[:12]),('tail_r',revcomp(OVERLAP)+revcomp(LEFT[-12:])),('tail_f',OVERLAP+RIGHT[:12]),('outer_r',revcomp(RIGHT[-12:]))):
        primers[identifier]=g.input(identifier,'invented '+identifier,K.OLIGO,state=dna(identifier,sequence),inferred=True,source=EXAMPLE)
    def amplify(identifier,template,forward,reverse):
        con=template.chemical_state
        geom=PCRGeometry.exact(con,PCRPrimer(forward.id,forward.state.top(),12,inferred=True),PCRPrimer(reverse.id,reverse.state.top(),12,inferred=True),statement='Exact binding in the invented example; minimum binding span is a modelling declaration, not an experimental recommendation.')
        product=pcr_product(geom,forward.state,reverse.state,name=identifier+' reference')
        out=g.product(identifier,'INFERRED — '+identifier+' reference',K.MOLECULAR,state=duplex(product),chemical_state=product)
        g.react(identifier+'_pcr',('INFERRED — illustrative outer-primer PCR; full original methods unreviewed' if identifier=='final' else 'INFERRED — amplify the declared toy reference'),S.PCR,[(template,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(pol,R.CATALYST)],[out],geometry=geom,relations=(M.PCR_PRODUCT,))
        return out
    af=amplify('fragment_a',a,primers['outer_f'],primers['tail_r'])
    bf=amplify('fragment_b',b,primers['tail_f'],primers['outer_r'])
    geom=OverlapExtension('fragment_a','fragment_b',af.chemical_state,bf.chemical_state,len(OVERLAP),'Exact toy end overlap; strand choice and copied reference are explicit. Denaturation allows this active strand pair; no pairing fraction or yield is assigned.')
    joined=geom.product('toy linear fusion')
    fused=g.product('fused','INFERRED — linear fused reference',K.MOLECULAR,state=duplex(joined),chemical_state=joined)
    # The active-strand view is the pre-extension molecular state, not a separate reaction.
    active=g.product('active','INFERRED — selected active-strand overlap pair',K.MOLECULAR,state=Drawing(tuple(geom.rows())),chemical_state=geom)
    g.react('pair','INFERRED — denature/reanneal to the illustrated productive pair',S.HYBRIDIZATION,[(af,R.TEMPLATE),(bf,R.STRUCTURAL)],[active],reason='Selected productive pair only; alternative strand pairings and residual primers remain outside this copying example.')
    g.react('extend','INFERRED — copy both flanks across the single overlap',S.EXTENSION,[(active,R.TEMPLATE),(pol,R.CATALYST)],[fused],geometry=geom,relations=(M.OVERLAP_EXTENSION,))
    final=amplify('final',fused,primers['outer_f'],primers['outer_r'])
    g.constructed(final,'Linear illustrative DNA reference; circular cloning, NGS preparation and experimental success are outside scope.')
    g.declare_reviewed_scope(('a','b','polymerase','outer_f','tail_r','tail_f','outer_r','fragment_a_pcr','fragment_b_pcr','pair','extend','final_pcr','final'),SOURCE+' Exact DNA, primer tails and binding sites are invented; the final outer-primer PCR is an explicitly chosen illustrative continuation, not established from the defining abstract.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_soe_pcr.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1016/0378-1119(89)90359-4">Horton et al.1989</a>. INFERRED — invented templates, primers and exact overlap. PCR reference products omit unassigned terminal heterogeneity, copying errors and competing species.'),
        WorkflowSection('Two tailed PCRs → overlap extension → outer-primer PCR'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory reference-copying claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
