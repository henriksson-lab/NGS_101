"""Berrow2007 In-Fusion Dry-Down: neutral reference, not a current-kit recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from pcr import PCRGeometry,PCRPrimer
from assembly_geometry import pcr_product,CircularOverlapReference
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='in_fusion.html'
TITLE='In-Fusion 2007 — illustrative Dry-Down implementation'
SOURCE='Berrow et al.2007 doi10.1093/nar/gkm047; published Methods In-Fusion cloning of PCR products and Discussion. The historical In-Fusion Dry-Down 96-well reagent is identified; current HD/Snap products and expression applications are outside scope.'
EXAMPLE='Invented short neutral DNA, primer tails and prepared recipient. No pOPIN map, source insert, functional vector, experimental performance or current-kit equivalence is claimed.'
LEFT='ACGTCAGTACGATCG'
RIGHT='TCGAGTACCTGACAT'
BODY='GATCTACGAGTCTCGATGACTAGC'
RECIPIENT='ATCGTAGCTCAGTCGATCAG'

def dna(name,seq):return Construct((Segment(name,seq,inferred=True),),name)
def limitation(text):return Drawing((Row(chunks=[(text,None,False)]),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    template=g.input('template','invented neutral insert template',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(dna('neutral template',BODY))),chemical_state=dna('neutral template',BODY),inferred=True,source=EXAMPLE)
    recipient_con=Construct((Segment('right insert homology',RIGHT,inferred=True),Segment('neutral recipient body',RECIPIENT,inferred=True),Segment('left insert homology',LEFT,inferred=True)),'invented prepared recipient reference')
    recipient=g.input('recipient','invented prepared linear recipient',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(recipient_con)),chemical_state=recipient_con,inferred=True,source=EXAMPLE+' Source recipient restriction digestion and purification precede this prepared-input boundary; no particular pOPIN cutting map is adopted.')
    primers=[]
    for identifier,sequence in (('forward',LEFT+BODY[:12]),('reverse',revcomp(RIGHT)+revcomp(BODY[-12:]))):
        primers.append(g.input(identifier,'invented '+identifier+' primer',K.OLIGO,state=dna(identifier,sequence),inferred=True,source=EXAMPLE))
    polymerase=g.input('polymerase','source KOD Hi-Fi polymerase',K.CATALYST)
    geom=PCRGeometry.exact(template.chemical_state,PCRPrimer('forward',primers[0].state.top(),12,inferred=True),PCRPrimer('reverse',primers[1].state.top(),12,inferred=True),statement='Exact illustrative primer sites and tails; the selected binding-span length is a model declaration, not a source recommendation.')
    insert_con=pcr_product(geom,primers[0].state,primers[1].state,name='invented tailed PCR reference')
    amplified=g.product('amplified','INFERRED — tailed insert PCR reference',K.MOLECULAR,state=Scene.duplex(tuple(insert_con)),chemical_state=insert_con)
    g.react('pcr','INFERRED — copy the neutral insert with 15-base homology tails',S.PCR,[(template,R.TEMPLATE),*( (p,R.STRUCTURAL) for p in primers),(polymerase,R.CATALYST)],[amplified],geometry=geom,relations=(M.PCR_PRODUCT,))
    purified=g.product('purified','INFERRED — purified insert reference',K.MOLECULAR,state=Scene.duplex(tuple(insert_con)),chemical_state=insert_con)
    waste=g.product('purification_waste','unexpanded PCR-purification waste',K.WASTE,state=limitation('Source PCR cleanup; exact primer, template and mixture partition is not assigned.'))
    g.react('purify','Purify the PCR insert before the In-Fusion reaction',S.CLEANUP,[(amplified,R.TEMPLATE)],[purified,waste],reason='Published PCR purification step; no quantitative removal or purity is assigned.')
    g.discard(waste,'Source purification waste; composition and removal efficiency unassigned.')
    reagent=g.input('infusion','historical In-Fusion Dry-Down 96-well reagent',K.CATALYST)
    host=g.input('host','source bacterial transformation/recovery host',K.BIOLOGICAL)
    reference=CircularOverlapReference((('purified',insert_con,len(LEFT),len(RIGHT)),('recipient',recipient_con,len(RIGHT),len(LEFT))),'Only the two declared 15-base reference overlaps are compared. The exact reagent-treated molecule and intracellular closure pathway are not assigned.')
    con=reference.reference('intended neutral circular reference')
    final=g.product('final','INFERRED — intended recovered circular reference',K.MOLECULAR,state=closed_circle('intended toy reference',tuple(con)),chemical_state=con)
    residual=g.product('residual','unresolved residual substrates and recovery alternatives',K.MOLECULAR,state=limitation('Unreacted DNA, alternative assemblies, reagent intermediates and host-associated fates remain unexpanded.'))
    g.react('infusion_and_recovery','INFERRED — historical Dry-Down treatment → bacterial recovery',S.MAGIC,[(purified,R.TEMPLATE),(recipient,R.STRUCTURAL),(reagent,R.CATALYST),(host,R.TEMPLATE)],[final,residual],geometry=reference,reason='Ordered source boundary: the prepared insert and recipient enter the historical In-Fusion reagent, then the diluted mixture is transformed into bacteria. No reagent-specific exonuclease/strand-exchange intermediate, covalent closure state or host-repair pathway is invented. A matching intended reference is not reaction validation.')
    g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'No conversion, clone recovery, background fraction or purity claimed.'))
    g.constructed(final,'Conditional intended constructed-DNA reference; biological function, expression, replication and sequencing are outside this neutral example.')
    g.declare_reviewed_scope(('template','recipient','forward','reverse','polymerase','pcr','purify','infusion','host','infusion_and_recovery','final'),SOURCE+' Prepared linear recipient boundary and exact example bases are declared locally; 15-base homology is supported by the source Discussion. Physical reagent and host intermediates remain unresolved.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_in_fusion.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/gkm047">Berrow et al.2007</a>. Historical In-Fusion Dry-Down implementation. INFERRED — neutral DNA, exact primer sites and intended reference. Reagent intermediates and bacterial repair remain unresolved; this is not a modern-kit recipe.'),
        WorkflowSection('Tailed PCR → purification → historical In-Fusion/recovery boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory declared reference claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkm047_PMC1874605.xml.txt','Published Methods and Discussion, historical implementation')
    if have(c,source,label='historical reagent identity and overlap scope'):
        text=source.path.read_text();c('source identifies Dry-Down and arbitrary 15bp homology','Dry-Down 96-well plate' in text and 'any 15' in text and 'homology region' in text,True)
    c.report()
