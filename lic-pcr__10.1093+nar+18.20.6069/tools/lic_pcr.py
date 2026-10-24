"""LIC-PCR: explicit invented nucleotide-boundary example, no experimental claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from pcr import PCRGeometry,PCRPrimer
from assembly_geometry import pcr_product,NucleotideBoundary
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='lic_pcr.html'
TITLE='LIC-PCR — illustrative nucleotide-controlled end processing'
SOURCE='Aslanidis and de Jong1990 doi10.1093/nar/18.20.6069; defining PubMed abstract PMID2235490. Primer-tail nucleotide exclusions, T4 DNA polymerase with dGTP/dCTP and transformation without in-vitro ligase; full original deferred.'
EXAMPLE='Invented short DNA and tailed primers, not used in the cited experiments. The chosen retention boundary is a composition example, not a prediction of enzymatic endpoints, cloning efficiency or a functioning vector.'
LEFT='AGTATGATAGTA'
RIGHT='ACTATCACTATA'
INSERT='CTAGTCAGATCGTACGATGCATAG'
RECIPIENT='GATCGTACAGCTAGTCGATGATCC'

def dna(name,seq):
    return Construct((Segment(name,seq,inferred=True),),name)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    pcr_pol=g.input('pcr_pol','source PCR DNA polymerase; detailed identity deferred',K.CATALYST)
    t4=g.input('t4','source T4 DNA polymerase exonuclease treatment',K.CATALYST)
    prepared=[];fragments=[]
    for identifier,core,left,right,nt,base in (('insert',INSERT,LEFT,RIGHT,'dGTP','G'),('recipient',RECIPIENT,RIGHT,LEFT,'dCTP','C')):
        con=dna(identifier+' core',core)
        parent=g.input(identifier+'_template','invented '+identifier+' template',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(con)),chemical_state=con,inferred=True,source=EXAMPLE)
        f=dna(identifier+' forward',left+core[:12]);r=dna(identifier+' reverse',revcomp(right)+revcomp(core[-12:]))
        fm=g.input(identifier+'_f','invented tailed primer',K.OLIGO,state=f,inferred=True,source=EXAMPLE)
        rm=g.input(identifier+'_r','invented tailed primer',K.OLIGO,state=r,inferred=True,source=EXAMPLE)
        pg=PCRGeometry.exact(con,PCRPrimer(f.name,f.top(),12,inferred=True),PCRPrimer(r.name,r.top(),12,inferred=True),statement='Exact invented PCR geometry; primer spans are model choices, not laboratory recommendations.')
        amplified=pcr_product(pg,f,r,name=identifier+' reference')
        amp=g.product(identifier+'_pcr_product','INFERRED — '+identifier+' tailed PCR reference',K.MOLECULAR,state=Scene.duplex(tuple(amplified)),chemical_state=amplified)
        g.react(identifier+'_pcr','INFERRED — generate toy complementary tails by PCR',S.PCR,[(parent,R.TEMPLATE),(fm,R.STRUCTURAL),(rm,R.STRUCTURAL),(pcr_pol,R.CATALYST)],[amp],geometry=pg,relations=(M.PCR_PRODUCT,))
        fragment=CohesiveFragment(identifier,left,(Segment(identifier+' retained core',core,inferred=True),),right,'5-prime','5-prime',inferred=True)
        geometry=NucleotideBoundary(amp.id,identifier+'_ready',amplified.top(),revcomp(amplified.top()),(0,len(amplified)-len(right)),(0,len(amplified)-len(left)),fragment,'Exact declared toy boundary: removed 3-prime suffixes exclude the supplied nucleotide base; first retained terminal base matches. This does not predict the actual exonuclease/polymerase equilibrium.',base)
        nucleotide=g.input(identifier+'_nucleotide',nt+' for the source '+identifier+' branch',K.CONDITION)
        ready=g.product(identifier+'_ready','INFERRED — processed '+identifier+' reference',K.MOLECULAR,state=fragment,chemical_state=fragment)
        removed=g.product(identifier+'_removed','removed parent reference intervals; chemical products unexpanded',K.MOLECULAR,state=Drawing((Row(chunks=[('Removed reference bases; individual degradation products and end-group chemistry are not drawn.',None,False)]),)))
        g.react(identifier+'_process','INFERRED — nucleotide-qualified 3-prime processing boundary',S.DIGESTION,[(amp,R.TEMPLATE),(t4,R.CATALYST),(nucleotide,R.CONDITION)],[ready,removed],geometry=geometry,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Removed-base chemical species and physical cleanup are unexpanded.'))
        prepared.append(ready);fragments.append(fragment)
    assembly=CircularCohesiveAssembly(tuple((m.id,f) for m,f in zip(prepared,fragments)),'Invented complementary tails and explicit circular order; no in-vitro ligase closure.')
    product=assembly.product('intended toy reference')
    joint=g.product('joint','INFERRED — unsealed toy assembly',K.MOLECULAR,state=Drawing(tuple(assembly.rows())),chemical_state=product)
    g.react('anneal','INFERRED — anneal without in-vitro ligase',S.HYBRIDIZATION,[(prepared[0],R.TEMPLATE),(prepared[1],R.STRUCTURAL)],[joint],geometry=assembly,relations=(M.COHESIVE_ASSEMBLY,))
    host=g.input('host','source bacterial transformation/recovery boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended repaired circular reference',K.MOLECULAR,state=closed_circle('toy circular reference',tuple(product)),chemical_state=product)
    g.react('recover','INFERRED — host repair boundary',S.MAGIC,[(joint,R.TEMPLATE),(host,R.TEMPLATE)],[final],reason='Intracellular repair is not reconstructed. The toy body has no replicon, so no viable clone or experimental recovery is claimed.')
    g.constructed(final,'Illustrative DNA reference, not a verified clone or NGS library.')
    g.declare_reviewed_scope(('pcr_pol','t4','insert_template','insert_f','insert_r','insert_pcr','insert_nucleotide','insert_process','recipient_template','recipient_f','recipient_r','recipient_pcr','recipient_nucleotide','recipient_process','anneal','host','recover','final'),SOURCE+' Exact templates, primers, nucleotide-boundary instances and body sequences are invented.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_lic_pcr.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/18.20.6069">Aslanidis and de Jong1990</a>. INFERRED — invented PCR substrates, primers and idealized retention boundaries. Supplied nucleotides constrain the illustrated composition; exact enzymatic stopping, terminal heterogeneity and yield are not established.'),
        WorkflowSection('Tailed PCRs → nucleotide-qualified processing → unsealed assembly → repair boundary'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory molecular claims in the declared example',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
