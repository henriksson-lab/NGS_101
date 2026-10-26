"""Source-qualified historical nanopore endpoint on a neutral example."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,Row
from assembly_geometry import DuplexPartition
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph,Disposition,Terminal,SequencingPlan,MaterialKind as K,InputRole as R,StepKind as S,MolecularRelation as M
from protocol_model import Drawing
from protocol_views import WorkflowSection,PlanLimitationSection,assemble_page
import seqprimers as sp
from page import head,info
PAGE_FILE='catch_nanopore.html'
TITLE='CATCH nanopore — illustrative native-DNA branch (2018)'
SOURCE='Gabrieli et al.2018 doi10.1093/nar/gky411; Methods Cas9 digestion, Separation by PFGE and Nanopore sequencing; Figure1. Native unamplified branch only.'
EXAMPLE='Invented short neutral reference window; not a complete source chromosome, published genomic target, guide design, functional construct or validated library. Flank intervals account for reference bases without predicting detached chromosome topology.'

def dna(name,bases):return Construct((Segment(name,bases,inferred=True),),name)
def duplex(con):return Scene.duplex(tuple(con),label=con.name)
def limitation(text):return Drawing((Row(chunks=[(text,None,False)]),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='illustrative_example',example_description=EXAMPLE)
    return build_graph(g)

def build_graph(g):
    target=dna('neutral target','ACGTCAGTACGATCGTGACTACGAGATCTATCGATGACTAGCATCGTACGATGCTAGTCA')
    parent=dna('neutral parent','GCTAGATC'+target.top()+'CTAGTCAG')
    source=g.input('parent','invented reference in gel context',K.MOLECULAR,primary=True,state=duplex(parent),inferred=True,source=EXAMPLE)
    guides=g.input('guides','source guide pair; example sequences and PAMs unassigned',K.OLIGO)
    nuclease=g.input('nuclease','source Cas9 nuclease',K.CATALYST)
    geom=DuplexPartition('parent',parent,(('left',0,8),('target',8,8+len(target)),('right',8+len(target),len(parent))),'Chosen reference boundaries; physical cut offsets, terminal groups and guide recognition remain unassigned.')
    cuts=[g.product(i,'INFERRED — '+i+' reference',K.MOLECULAR,state=duplex(con)) for i,con in geom.products()]
    g.react('cut','INFERRED — source delimitation; chosen reference boundaries',S.DIGESTION,[(source,R.TEMPLATE),(guides,R.STRUCTURAL),(nuclease,R.CATALYST)],cuts,geometry=geom,relations=(M.DUPLEX_PARTITION,))
    proteinase=g.input('proteinase','source post-cleavage proteinase K',K.CATALYST)
    treated=g.product('treated','INFERRED — target after deproteinization',K.MOLECULAR,state=duplex(target),chemical_state=target)
    treated_flanks=[g.product(i+'_treated','INFERRED — '+i+' after deproteinization',K.MOLECULAR,state=duplex(con),chemical_state=con) for i,con in geom.products() if i!='target']
    other_protein=g.product('protein_products','unexpanded protein/reagent products',K.MOLECULAR,state=limitation('Protein digestion products remain unexpanded.'))
    g.react('deproteinize','Remove bound/excess protein from the whole digest before separation',S.DEPROTEINIZATION,[(m,R.TEMPLATE) for m in cuts]+[(proteinase,R.CATALYST)],[treated,*treated_flanks,other_protein])
    g.mark_terminal(Terminal(other_protein.id,Disposition.UNRESOLVED,'Protein/reagent fates remain unexpanded.'))
    selection=g.input('pfge','source PFGE and target-band excision',K.SELECTION)
    band=g.product('band','INFERRED — selected reference',K.MOLECULAR,state=duplex(target),chemical_state=target)
    background=g.product('background','non-selected gel material',K.WASTE)
    g.react('isolate','PFGE and excise target-band branch',S.SELECTION,[(m,R.TEMPLATE) for m in (treated,*treated_flanks)]+[(selection,R.SELECTION)],[band,background],reason='Reference tracking does not predict migration, co-migration, purity or recovery.')
    g.discard(background,'Non-excised gel is physically left behind; complete background removal is not established.')
    recovery=g.input('recovery','source agarase and precipitation recovery',K.SELECTION)
    native=g.product('native','INFERRED — recovered native reference',K.MOLECULAR,state=duplex(target),chemical_state=target)
    lost=g.product('loss','unresolved recovery loss',K.MOLECULAR,state=limitation('Recovery losses and co-isolated DNA remain unassigned.'))
    g.react('recover','Recover DNA from the selected band',S.EXTRACTION,[(band,R.TEMPLATE),(recovery,R.SELECTION)],[native,lost])
    g.mark_terminal(Terminal(lost.id,Disposition.UNRESOLVED,'No yield or purity is inferred.'))
    return finish_graph(g,native)

def finish_graph(g,native):
    kit=g.input('kit','source low input expansion pack plus SQK-MAP007 R9',K.ADAPTER)
    unknown=AdapterLigationGeometry.unavailable('The article identifies SQK-MAP007 R9 plus a low input expansion pack. The applicable historical kit protocol and adapter-end geometry have not been reviewed; modern ONT chemistry is not substituted.')
    library=g.product('library','historical-kit library boundary',K.LIBRARY,state=limitation('Final adapter polymers, end geometry and molecule classes remain unassigned. The neutral DNA above is a substrate reference, not a resolved final-library drawing.'),inferred=True,source=SOURCE)
    other=g.product('kit_other','unresolved kit products and residual substrates',K.MOLECULAR,state=limitation('Adapter excess, losses and alternative products are unexpanded.'))
    g.react('prepare','Historical kit preparation boundary; source branch is unamplified',S.MAGIC,[(native,R.TEMPLATE),(kit,R.STRUCTURAL)],[library,other],geometry=unknown,reason='Kit identity and sequencing endpoint are source-established. Historical kit chemistry remains unresolved; no end-prep, ligase, motor or hairpin geometry is inferred.')
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Historical-kit byproducts and complete conversion remain unassigned.'))
    g.sequence(library,SequencingPlan(limitation='Published native-DNA branch: low input expansion pack plus SQK-MAP007 R9, sequenced on R9.4 flow cells. Adapter geometry, strand acquisition and per-molecule read length remain unassigned. Human MDA/T7-endonuclease/SQK-LSK308 R9.5 and companion Illumina branches are outside scope.',acquired_roles=('Sequencing output',),read_uses={'Sequencing output':'include'}))
    g.declare_reviewed_scope(('parent','guides','nuclease','cut','proteinase','deproteinize','pfge','isolate','recovery','recover','kit','prepare','library'),SOURCE+' Example DNA and boundaries are invented; historical kit protocol is deferred, not claimed absent.')
    return g.require_publishable()

def render_page():
    graph=protocol_graph()  # require_publishable precedes the completion marker.
    plan=graph.sequencing['library']
    return assemble_page(graph,(head(TITLE),'<div class="wrap" data-protocol-complete="true"><h1>'+TITLE+'</h1>','<p class="research-notes"><a href="01_catch_nanopore.html">Research notes</a></p>',info('<a href="https://doi.org/10.1093/nar/gky411">Gabrieli et al.2018</a>. INFERRED — neutral DNA and reference boundaries. Historical kit internals remain unresolved. The amplified human-target branch uses a different preparation.'),WorkflowSection('Delimitation → band isolation → recovery → historical kit → nanopore'),PlanLimitationSection('library'),sp.read_spec_table([sp.read_spec(role,None,'Native target DNA; historical-kit strand classes and adapter geometry remain unassigned.',plan.read_uses[role]) for role in plan.acquired_roles],'Source native-DNA branch, R9.4 nanopore readout; no synthesis sequencing primer or fixed cycle count.'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory reference claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
