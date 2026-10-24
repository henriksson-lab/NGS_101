"""BASIC linker workflow with explicitly customized, invented part ends."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from plasmid import BsaI,Plasmid,digest
from assembly_geometry import StrandRetention,LinearFivePrimeAssembly
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='basic.html'
TITLE='BASIC — illustrative customized linker assembly'
SOURCE='Storch et al.2015 doi10.1021/sb500356d; published publisher text/Figure1 caption: BsaI/linker ligation, excess-linker purification, annealing and host recovery. Main text permits adapting linker-ligation ends to other TypeIIS parts. Published iP/iS and actual neutral-linker tables remain unassigned; source SI inaccessible.'
EXAMPLE='Invented short DNA bodies, customized part-end/linker sequences and prepared end groups. These are not the published BASIC iP/iS standard, source oligos, an experimental instance or a working vector. The example illustrates the source linker-assembly mechanism and claims no standard-part compatibility.'
# Source lengths describe the general linker architecture; bases are invented.
J1='ACGTACGATCTAGTCAGATCG'
J2='TGCACTAGCATCGATGACTAG'
CUT_LEFT='AGTC'
CUT_RIGHT='CTGA'

def bp(name,seq):return Segment(name,seq,inferred=True)
def frag(name,left,core,right):return CohesiveFragment(name,left,tuple(core),right,'5-prime','5-prime',inferred=True)
def limit(text):return Drawing((Row(chunks=[(text,None,False)]),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    bsai=g.input('bsai','source BsaI',K.CATALYST);ligase=g.input('ligase','source T4 DNA ligase',K.CATALYST)
    beads=g.input('beads','source excess-linker purification',K.SELECTION)
    ready=[];adapted=[]
    for identifier,body,left_joint,right_joint in (('a','CGATGACCTGATAC',J1,J2),('b','TCAGGATCGATACG',J2,J1)):
        part=frag('toy part '+identifier,CUT_LEFT,(bp('body '+identifier,body),),CUT_RIGHT)
        sequence=BsaI.site+'A'+CUT_LEFT+body+CUT_RIGHT+'T'+revcomp(BsaI.site)
        parent=g.input(identifier+'_source','invented BsaI-flanked precursor',K.MOLECULAR,primary=True,state=Scene.duplex((bp('custom precursor',sequence),)),inferred=True,source=EXAMPLE)
        cuts=digest(Plasmid('toy precursor',sequence,circular=False),BsaI)
        if len(cuts)!=2 or (cuts[0].strand,cuts[1].strand)!=(1,-1):raise ValueError('custom precursor needs exactly the declared site pair')
        lo,hi=cuts
        geom=StrandRetention(parent.id,identifier+'_part',sequence,revcomp(sequence),(lo.top_cut,hi.top_cut),(len(sequence)-hi.bottom_cut,len(sequence)-lo.bottom_cut),part,'Invented cleavage cassette, not the canonical BASIC iP/iS; retained strands follow imported BsaI cut coordinates.')
        cut=g.product(identifier+'_part','INFERRED — released customized part',K.MOLECULAR,state=part,chemical_state=part)
        remnants=g.product(identifier+'_remnants','cleaved precursor remainder class',K.MOLECULAR,state=limit('Removed reference intervals; exact remainder duplexes and handling unexpanded.'))
        g.react(identifier+'_cut','INFERRED — BsaI release within linker-attachment stage',S.DIGESTION,[(parent,R.TEMPLATE),(bsai,R.CATALYST)],[cut,remnants],geometry=geom,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(remnants.id,Disposition.UNRESOLVED,'Precursor remainder fate is unexpanded; no complete clearance or ligation yield claimed.'))
        prefix=frag('invented prefix linker',left_joint,(bp('prefix paired core','GACCTAGTCGAT'),),CUT_LEFT)
        suffix=frag('invented suffix linker',CUT_RIGHT,(bp('suffix paired core','CTAGACGTCAGT'),),right_joint)
        pm=g.input(identifier+'_prefix','invented paired prefix-linker oligos',K.OLIGO,state=prefix,chemical_state=prefix,inferred=True,source=EXAMPLE)
        sm=g.input(identifier+'_suffix','invented paired suffix-linker oligos',K.OLIGO,state=suffix,chemical_state=suffix,inferred=True,source=EXAMPLE)
        geometry=LinearFivePrimeAssembly(((pm.id,prefix),(cut.id,part),(sm.id,suffix)),True,True,'Invented prepared joining ends explicitly have 5-prime phosphates and 3-prime hydroxyls. Canonical source oligo modifications and iP/iS sequences are not assigned.')
        product=geometry.product('linker-adapted '+identifier)
        attached=g.product(identifier+'_attached','INFERRED — linker-adapted toy part',K.MOLECULAR,state=product,chemical_state=product)
        free=g.product(identifier+'_free','unligated linker class in reaction mixture',K.MOLECULAR,state=limit('Unligated linker population; abundance and sequences of mixture species unassigned.'))
        g.react(identifier+'_attach','INFERRED — attach customized linkers in source restriction–ligation stage',S.LIGATION,[(cut,R.TEMPLATE),(pm,R.STRUCTURAL),(sm,R.STRUCTURAL),(ligase,R.CATALYST)],[attached,free],geometry=geometry,relations=(M.LINEAR_COHESIVE_ASSEMBLY,))
        purified=g.product(identifier+'_purified','INFERRED — selected linker-adapted part class',K.MOLECULAR,state=product,chemical_state=product)
        removed=g.product(identifier+'_removed','separated excess-linker class',K.WASTE)
        g.react(identifier+'_cleanup','Separate excess unligated linkers',S.CLEANUP,[(attached,R.TEMPLATE),(free,R.TEMPLATE),(beads,R.SELECTION)],[purified,removed])
        g.discard(removed,'Source purification removes excess oligos. A selected class is illustrated; quantitative recovery or perfect purity is not asserted.')
        ready.append(purified);adapted.append(product)
    assembly=CircularCohesiveAssembly(tuple((m.id,f) for m,f in zip(ready,adapted)),'Customized long linker ends determine the declared circular order; both linker cores are retained, each long overlap once.')
    con=assembly.product('intended customized reference')
    joint=g.product('joint','INFERRED — unsealed linker-mediated assembly',K.MOLECULAR,state=Drawing(tuple(assembly.rows())),chemical_state=con)
    g.react('anneal','INFERRED — anneal adapted parts without added ligase',S.HYBRIDIZATION,[(ready[0],R.TEMPLATE),(ready[1],R.STRUCTURAL)],[joint],geometry=assembly,relations=(M.COHESIVE_ASSEMBLY,))
    host=g.input('host','source host recovery boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended host-repaired reference',K.MOLECULAR,state=closed_circle('intended toy reference',tuple(con)),chemical_state=con)
    g.react('recover','INFERRED — host repair boundary',S.MAGIC,[(joint,R.TEMPLATE),(host,R.TEMPLATE)],[final],reason='Intracellular repair and actual cloning are outside scope; the short toy has no replicon.')
    g.constructed(final,'Illustrative customized DNA reference, not standard BASIC compatibility, a verified clone or an NGS library.')
    g.declare_reviewed_scope(('bsai','ligase','beads','a_source','a_cut','a_prefix','a_suffix','a_attach','a_cleanup','b_source','b_cut','b_prefix','b_suffix','b_attach','b_cleanup','anneal','host','recover','final'),SOURCE+' Customized bases, cleavage cassettes and end-preparation preconditions are invented, not source transcription.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_basic.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1021/sb500356d">Storch et al.2015</a>. INFERRED — customized example part ends, linkers and supplied end groups. Canonical BASIC iP/iS and source oligos are unassigned. Molecular panels split relations within the source workflow, not separate experimental incubations.'),
        WorkflowSection('Release and linker attachment → excess-linker removal → annealing → repair boundary'),'</div>'))

def selftest():
    from checks import Check
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory customized-example relation',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[]);c.report()
