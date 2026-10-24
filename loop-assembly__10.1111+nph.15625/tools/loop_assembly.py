"""Published Loop alternation; invented receiver/fusion grammar, no source maps."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from plasmid import BsaI,Enzyme,Plasmid,digest
from assembly_geometry import StrandRetention,dna_window
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='loop_assembly.html'
TITLE='Loop assembly — illustrative two-level recursion'
SOURCE='Pollak et al.2019 doi10.1111/nph.15625; published HTML Introduction, Loop assembly Results and Figures1–2 captions. SapI recognition/cut offsets independently from NEB recognition-specificity table. Superseded preprint supplements excluded.'
EXAMPLE='Invented neutral DNA bodies, receiver maps and fusion grammar. These are not the published pOdd/pEven sequences, PhytoBrick standard parts or an experimental instance; no functional vector or cloning performance is claimed.'
SAPI=Enzyme('SapI','GCTCTTC',1,4)
ODD_ENDS=('AGTC','CTGA','TCGT','GACT','ACGA')
EVEN_ENDS=('ATC','GAC','TCA','CGT','GTA')
BODIES=('CGATGACC','TCAGGATC','ATCGTAGC','GATCAGTC')

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
    bsai=g.input('bsai','source odd-level BsaI',K.CATALYST)
    sapi=g.input('sapi','source even-level SapI',K.CATALYST)
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
        geometry=CircularCohesiveAssembly(tuple((m.id,f) for m,f in parts),'Declared toy circular order; future alternate-enzyme sites remain inside the retained receiver body.')
        con=rotate_at(geometry.product(identifier+' circular reference'),next_site)
        out=g.product(identifier,'INFERRED — '+identifier+' circular reference',K.MOLECULAR,state=closed_circle(identifier+' toy reference',tuple(con)),chemical_state=Scene.duplex(tuple(con)))
        g.react(identifier+'_ligate','INFERRED — ligate the declared '+identifier+' order',S.LIGATION,[(parts[0][0],R.TEMPLATE)]+[(m,R.STRUCTURAL) for m,_ in parts[1:]]+[(ligase,R.CATALYST)],[out],geometry=geometry,relations=(M.COHESIVE_ASSEMBLY,))
        return out
    # The four input classes can contribute to parallel example assemblies; this
    # is not a declaration of physical aliquots, amplification or copy number.
    l0=[]
    for i,body in enumerate(BODIES):
        con=precursor(BsaI,ODD_ENDS[i],(bp('neutral part '+str(i+1),body),),ODD_ENDS[i+1],'L0 toy')
        l0.append(release('l0_'+str(i+1),source('l0_source_'+str(i+1),con),BsaI,bsai))
    l1=[]
    for i in range(4):
        core=(bp('future SapI right end',EVEN_ENDS[i+1]),bp('future right spacer','T'),bp('future SapI reverse',revcomp(SAPI.site)),bp('receiver body','TCGATAGC'),bp('future SapI forward',SAPI.site),bp('future left spacer','A'),bp('future SapI left end',EVEN_ENDS[i]))
        receiver=source('odd_source_'+str(i+1),precursor(BsaI,ODD_ENDS[-1],core,ODD_ENDS[0],'odd receiver toy'))
        receiver_part=release('odd_receiver_'+str(i+1),receiver,BsaI,bsai)
        l1.append(ligate('l1_'+str(i+1),[*l0,receiver_part],'future SapI forward'))
    l1_parts=[release('l1_donor_'+str(i+1),m,SAPI,sapi) for i,m in enumerate(l1)]
    # Retaining the next BsaI pair is the physical return path to the odd set.
    core=(bp('future BsaI right end','TGAC'),bp('future right spacer','T'),bp('future BsaI reverse',revcomp(BsaI.site)),bp('even receiver body','GATCATGC'),bp('future BsaI forward',BsaI.site),bp('future left spacer','A'),bp('future BsaI left end','AGCT'))
    even=source('even_source',precursor(SAPI,EVEN_ENDS[-1],core,EVEN_ENDS[0],'even receiver toy'))
    receiver_part=release('even_receiver',even,SAPI,sapi)
    final=ligate('l2_final',[*l1_parts,receiver_part],'future BsaI forward')
    reference=Construct(tuple(final.chemical_state.strands['top'].segs))
    future=digest(Plasmid('final reference',reference.top(),circular=True),BsaI)
    if len(future)!=2 or any(len(c.overhang)!=4 for c in future):raise ValueError('return path must retain the intended future BsaI pair')
    g.constructed(final,'Conditional two-level toy reference; retained BsaI sites permit a further declared-map release. No actual level3 assembly, functional vector or source-standard compatibility claimed.')
    g.declare_reviewed_scope(('bsai','sapi','ligase','l0_source_1','l0_source_2','l0_source_3','l0_source_4','l0_1_cut','l0_2_cut','l0_3_cut','l0_4_cut','odd_source_1','odd_source_2','odd_source_3','odd_source_4','odd_receiver_1_cut','odd_receiver_2_cut','odd_receiver_3_cut','odd_receiver_4_cut','l1_1_ligate','l1_2_ligate','l1_3_ligate','l1_4_ligate','l1_donor_1_cut','l1_donor_2_cut','l1_donor_3_cut','l1_donor_4_cut','even_source','even_receiver_cut','l2_final_ligate','l2_final'),SOURCE+' Exact neutral maps and fusion strings are invented; receptor positions and biological bodies are not transcribed from source maps.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_loop_assembly.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1111/nph.15625">Published Pollak et al.2019</a>. INFERRED — invented maps and fusion grammar. BsaI/SapI identities and alternating architecture follow the published method. <a href="https://www.neb.com/en-us/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities">NEB</a> independently specifies SapI GCTCTTC(1/4).'),
        info('Published wording discrepancy: Results and Figure1a assign BsaI to odd receivers and SapI to even receivers. The publisher HTML Figure2b caption calls its BsaI receiver both odd and even. This model follows Results/Figure1a; the intended caption correction is unconfirmed. Superseded preprint maps are not used to resolve it.'),
        WorkflowSection('BsaI odd-level assembly → SapI even-level assembly → retained BsaI return sites'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    from chemistry_validation import check_graph,Status
    c=Check();c('no contradictory example strand/assembly claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'alphabetized-list-of-recognition-specificities.html.txt','NEB primary specificity table; browser verified, local retrieval403')
    if have(c,source,label='local supplier recognition/cut transcription'):
        text=''.join(source.path.read_text().split())
        c('source SapI motif and declared cut offsets',SAPI.site+'('+str(SAPI.top_cut)+'/'+str(SAPI.bottom_cut)+')' in text,True)
    c.report()
