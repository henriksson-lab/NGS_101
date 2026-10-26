"""Original GoldenBraid2011: neutral two-loop reference, no GB2.0 grammar."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from plasmid import BsaI,BsmBI,Plasmid,digest
from assembly_geometry import StrandRetention,dna_window
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='goldenbraid.html'
TITLE='GoldenBraid 2011 — illustrative binary double loop'
SOURCE='Sarrion-Perdigones et al.2011 doi10.1371/journal.pone.0021622; published Results Part standardization/The Double loop Design, Figures1–3 captions, Methods Cloning and assembly of modular pieces. Original BsaI/BsmBI four-receiver system, not GoldenBraid2.0. Supplier recognition table/product identity consulted only for cleavage geometry, not reaction formulations.'
EXAMPLE='Invented short neutral bodies, receiver references and four-base fusion strings. The example illustrates original2011 topology, not published pDGB maps, GB2.0 compatibility, coding devices, a functional vector or experimental performance.'
# The labels here identify invented strings; none is a source-standard fusion.
ONE='AGTC';TWO='GACT';THREE='TCGT'
A='ATCG';B='CAGT';C='GTAC'
INNER=('CTGA','ACGA')
BODIES=('CGATGACC','TCAGGATC','ATCGTAGC')

def bp(name,seq):return Segment(name,seq,inferred=True)
def limitation(text):return Drawing((Row(chunks=[(text,None,False)]),))
def precursor(enzyme,left,core,right,name):
    return Construct((bp('cut forward',enzyme.site),bp('left spacer','A'),bp('left end',left),*core,bp('right end',right),bp('right spacer','T'),bp('cut reverse',revcomp(enzyme.site)),bp('unretained body','CATAGTCA')),name)
def rotate_at(con,suffix):
    parts=tuple(con);indices=[i for i,s in enumerate(parts) if s.name.endswith(suffix)]
    if len(indices)!=1:raise ValueError('future site must have one declared segment identity')
    j=indices[0];return Construct(parts[j:]+parts[:j],con.name)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    bsai=g.input('bsai','source level-alpha BsaI',K.CATALYST)
    bsmbi=g.input('bsmbi','source level-omega BsmBI',K.CATALYST)
    ligase=g.input('ligase','source T4 DNA ligase',K.CATALYST)
    def source(identifier,con):
        return g.input(identifier,'invented '+identifier+' circular reference',K.MOLECULAR,primary=True,state=closed_circle('toy precursor',tuple(con)),chemical_state=Scene.duplex(tuple(con)),inferred=True,source=EXAMPLE)
    def release(identifier,parent,enzyme,catalyst):
        scene=parent.chemical_state
        con=Construct(tuple(scene.strands['top'].segs),parent.id+' opened reference')
        sequence=con.top();cuts=digest(Plasmid('toy circular map',sequence,circular=True),enzyme)
        if len(cuts)!=2 or (cuts[0].strand,cuts[1].strand)!=(1,-1):raise ValueError('toy reference origin must expose exactly the declared inward-cutting pair')
        lo,hi=cuts
        left=sequence[lo.top_cut:lo.bottom_cut];right=sequence[hi.top_cut:hi.bottom_cut]
        fragment=CohesiveFragment(identifier,left,dna_window(con,lo.bottom_cut,hi.top_cut,prefix='kept '),right,'5-prime','5-prime',inferred=True)
        geometry=StrandRetention(parent.id,identifier,sequence,revcomp(sequence),(lo.top_cut,hi.top_cut),(len(con)-hi.bottom_cut,len(con)-lo.bottom_cut),fragment,'Invented circular reference is opened outside retained cuts. Exact strand intervals are checked; enzymatic efficiency and remainder molecules are not predicted.')
        out=g.product(identifier,'INFERRED — '+identifier+' released reference',K.MOLECULAR,state=fragment,chemical_state=fragment)
        remainder=g.product(identifier+'_remainder','unretained precursor-reference class',K.MOLECULAR,state=limitation('Unretained parent intervals and mixture fates remain unexpanded.'))
        g.react(identifier+'_cut','INFERRED — '+enzyme.name+' release on the declared toy map',S.DIGESTION,[(parent,R.TEMPLATE),(catalyst,R.CATALYST)],[out,remainder],geometry=geometry,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'No physical removal, complete digestion or quantitative product partition asserted.'))
        return out,fragment
    def ligate(identifier,parts,next_site):
        geometry=CircularCohesiveAssembly(tuple((m.id,f) for m,f in parts),'Declared neutral circular order. Recognition sites of the current enzyme are excluded while the alternate-enzyme pair remains in the retained receiver body.')
        con=rotate_at(geometry.product(identifier+' circular reference'),next_site)
        out=g.product(identifier,'INFERRED — '+identifier+' circular reference',K.MOLECULAR,state=closed_circle(identifier+' toy reference',tuple(con)),chemical_state=Scene.duplex(tuple(con)))
        g.react(identifier+'_ligate','INFERRED — ligate the declared '+identifier+' order',S.LIGATION,[(parts[0][0],R.TEMPLATE)]+[(m,R.STRUCTURAL) for m,_ in parts[1:]]+[(ligase,R.CATALYST)],[out],geometry=geometry,relations=(M.COHESIVE_ASSEMBLY,))
        return out
    def receiver(identifier,current,outer_left,outer_right,future,next_left,next_right):
        # Current-enzyme cut sites bracket the retained backbone. The alternate
        # pair faces inward across the assembled insert after circular ligation.
        core=(bp('future '+future.name+' right end',next_right),bp('future right spacer','T'),bp('future '+future.name+' reverse',revcomp(future.site)),bp('neutral receiver body','TCGATAGC'),bp('future '+future.name+' forward',future.site),bp('future left spacer','A'),bp('future '+future.name+' left end',next_left))
        parent=source(identifier,precursor(current,outer_right,core,outer_left,identifier+' prepared receiver'))
        return release(identifier+'_kept',parent,current,bsai if current is BsaI else bsmbi)
    # Four alpha outputs supply the two complementary alpha receiver classes
    # to each of two omega assemblies. These are reference classes, not a claim
    # about experimental aliquots, cloning yield or material copy number.
    alpha=[]
    for i in range(4):
        parts=[];ends=(ONE,*INNER,TWO)
        for j,body in enumerate(BODIES):
            identifier='part_'+str(i+1)+'_'+str(j+1)
            parent=source(identifier+'_source',precursor(BsaI,ends[j],(bp('neutral part '+str(j+1),body),),ends[j+1],identifier+' toy precursor'))
            parts.append(release(identifier,parent,BsaI,bsai))
        left,right=((A,C) if i%2==0 else (C,B))
        rec=receiver('alpha_receiver_'+str(i+1),BsaI,ONE,TWO,BsmBI,left,right)
        alpha.append(ligate('alpha_'+str(i+1),[*parts,rec],'future '+BsmBI.name+' forward'))
    donors=[release('alpha_donor_'+str(i+1),m,BsmBI,bsmbi) for i,m in enumerate(alpha)]
    omega=[]
    for i,(left,right) in enumerate(((ONE,THREE),(THREE,TWO))):
        rec=receiver('omega_receiver_'+str(i+1),BsmBI,A,B,BsaI,left,right)
        omega.append(ligate('omega_'+str(i+1),[donors[2*i],donors[2*i+1],rec],'future '+BsaI.name+' forward'))
    omega_donors=[release('omega_donor_'+str(i+1),m,BsaI,bsai) for i,m in enumerate(omega)]
    rec=receiver('return_alpha_receiver',BsaI,ONE,TWO,BsmBI,A,C)
    final=ligate('return_alpha',[*omega_donors,rec],'future '+BsmBI.name+' forward')
    reference=Construct(tuple(final.chemical_state.strands['top'].segs))
    future=digest(Plasmid('final toy circular reference',reference.top(),circular=True),BsmBI)
    if len(future)!=2 or any(len(c.overhang)!=4 for c in future):raise ValueError('binary return must retain exactly the intended alternate-enzyme pair')
    g.constructed(final,'Conditional ligated neutral DNA reference after alpha → omega → alpha. Host propagation, selection, functional sequence and further assembly are outside scope.')
    reviewed=('bsai', 'bsmbi', 'ligase', 'return_alpha', 'part_1_1_source', 'part_1_1_cut', 'part_1_2_source', 'part_1_2_cut', 'part_1_3_source', 'part_1_3_cut', 'alpha_receiver_1', 'alpha_receiver_1_kept_cut', 'alpha_1_ligate', 'alpha_donor_1_cut', 'part_2_1_source', 'part_2_1_cut', 'part_2_2_source', 'part_2_2_cut', 'part_2_3_source', 'part_2_3_cut', 'alpha_receiver_2', 'alpha_receiver_2_kept_cut', 'alpha_2_ligate', 'alpha_donor_2_cut', 'part_3_1_source', 'part_3_1_cut', 'part_3_2_source', 'part_3_2_cut', 'part_3_3_source', 'part_3_3_cut', 'alpha_receiver_3', 'alpha_receiver_3_kept_cut', 'alpha_3_ligate', 'alpha_donor_3_cut', 'part_4_1_source', 'part_4_1_cut', 'part_4_2_source', 'part_4_2_cut', 'part_4_3_source', 'part_4_3_cut', 'alpha_receiver_4', 'alpha_receiver_4_kept_cut', 'alpha_4_ligate', 'alpha_donor_4_cut', 'omega_receiver_1', 'omega_receiver_1_kept_cut', 'omega_1_ligate', 'omega_donor_1_cut', 'omega_receiver_2', 'omega_receiver_2_kept_cut', 'omega_2_ligate', 'omega_donor_2_cut', 'return_alpha_receiver', 'return_alpha_receiver_kept_cut', 'return_alpha_ligate')
    g.declare_reviewed_scope(reviewed,SOURCE+' All exact reference bodies, fusion strings and prepared receiver maps are invented. Current-enzyme sites are removed and the next enzyme sites retained by explicit strand-coordinate declarations; source functional assemblies, BbsI domestication and host recovery are outside this scope.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_goldenbraid.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1371/journal.pone.0021622">Original GoldenBraid2011</a>. INFERRED — neutral receiver maps and fusion strings. The source four-receiver alpha/omega topology uses BsaI/BsmBI; exact standard-map compatibility and biological performance are not claimed.'),
        WorkflowSection('Multipartite alpha assemblies → binary omega assemblies → binary alpha return'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory declared cut/assembly claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'journal.pone.0021622_PMC3131274.xml.txt','Original2011 Double loop and assembly methods')
    if have(c,source,label='original publication identity and four-receiver scope'):
        text=source.path.read_text();c('source preserves both original alpha and omega receiver pairs',all(s in text for s in ('pDGB_A12C','pDGB_C12B','pDGB_1AB3','pDGB_3AB2','T4 DNA ligase')),True)
    c.report()
