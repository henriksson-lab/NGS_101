"""Original CATCH cloning: neutral reference example, no functional cargo."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,Row,revcomp
from assembly_geometry import DuplexPartition,StrandRetention
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from protocol_graph import ProtocolGraph,Disposition,Terminal,MaterialKind as K,InputRole as R,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='catch_cloning.html'
TITLE='CATCH cloning — illustrative neutral capture'
SOURCE='Jiang et al.2015 doi10.1038/ncomms9101; Results CATCH and cloning, Methods DNA preparation and Vector preparation and cloning, Figure1 caption. Original cloning branch, not the 2016 protocol revision or nanopore adaptation.'
EXAMPLE='Invented short neutral reference window and prepared recipient, not complete source chromosomes, genomic targets, a functional vector or an experimental design. Flank intervals are reference accounting, not predictions of detached chromosome topology. Target boundaries, overlaps and exact resection endpoints are illustrative.'
LEFT='ACGTCAGTACGATCGTGACTACGAGATCTA'
RIGHT='TCGATGACTAGCATCGTACGATGCTAGTCA'
BODY='GATCGATACGTAGCTACGTA'
RECIPIENT='CATGACTAGTCAGCTAGCTA'

def dna(name,bases):return Construct((Segment(name,bases,inferred=True),),name)
def duplex(con):return Scene.duplex(tuple(con),label=con.name)
def limitation(text):return Drawing((Row(chunks=[(text,None,False)]),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    # All DNA is neutral; cleavage positions are declared reference coordinates.
    target=dna('neutral target reference',LEFT+BODY+RIGHT)
    parent=dna('neutral parent reference','GCTAGATC'+target.top()+'CTAGTCAG')
    substrate=g.input('parent','invented DNA reference in source-style gel context',K.MOLECULAR,primary=True,state=duplex(parent),inferred=True,source=EXAMPLE)
    guides=g.input('guides','source paired guides; sequences and PAM coordinates unassigned in this example',K.OLIGO,source=SOURCE)
    cas9=g.input('cas9','source Cas9 nuclease',K.CATALYST)
    geom=DuplexPartition('parent',parent,(('left_flank',0,8),('target',8,8+len(target)),('right_flank',8+len(target),len(parent))),'Chosen neutral reference boundaries; this is not a Cas9 cleavage-site or end-geometry prediction.')
    cuts=[g.product(i,'INFERRED — '+i+' reference',K.MOLECULAR,state=duplex(con)) for i,con in geom.products()]
    g.react('cut','INFERRED — source Cas9 delimitation; chosen reference boundaries',S.DIGESTION,[(substrate,R.TEMPLATE),(guides,R.STRUCTURAL),(cas9,R.CATALYST)],cuts,geometry=geom,relations=(M.DUPLEX_PARTITION,))
    gel=g.input('gel_recovery','source agarase and DNA recovery reagents',K.SELECTION)
    recovered=g.product('recovered','INFERRED — target reference after gel dissolution/recovery',K.MOLECULAR,state=duplex(target))
    recovery_other=g.product('recovery_other','unexpanded recovery background',K.MOLECULAR,state=limitation('Source recovery handles digested plug DNA without preparative target-band isolation. Background and losses are not quantified.'))
    g.react('recover','Dissolve gel and recover DNA; target reference follows the retained branch',S.EXTRACTION,[(m,R.TEMPLATE) for m in cuts]+[(gel,R.SELECTION)],[recovered,recovery_other])
    g.mark_terminal(Terminal(recovery_other.id,Disposition.UNRESOLVED,'Background and recovery losses remain unassigned; reference tracking does not imply target purity.'))
    recipient=dna('prepared neutral recipient',RIGHT+RECIPIENT+LEFT)
    rec=g.input('recipient','invented prepared recipient reference',K.MOLECULAR,primary=True,state=duplex(recipient),inferred=True,source=EXAMPLE+' Source PCR/DpnI preparation is outside this prepared-input scope.')
    exonuclease=g.input('exonuclease','source Gibson mixture T5 exonuclease',K.CATALYST)
    ready=[]
    for identifier,source,con,left,body,right in (('target',recovered,target,LEFT,BODY,RIGHT),('recipient',rec,recipient,RIGHT,RECIPIENT,LEFT)):
        fragment=CohesiveFragment(identifier,left,(Segment('neutral body',body,inferred=True),),right,'3-prime','3-prime',inferred=True)
        reference=con.top()
        retention=StrandRetention(source.id,identifier+'_ready',reference,revcomp(reference),(len(left),len(reference)),(len(right),len(reference)),fragment,'Exact resection endpoints chosen for the neutral example; enzyme stopping distributions and original cut-end offsets are not inferred.')
        out=g.product(identifier+'_ready','INFERRED — declared processed '+identifier,K.MOLECULAR,state=fragment)
        removed=g.product(identifier+'_removed','unexpanded removed '+identifier+' intervals',K.MOLECULAR,state=limitation('Removed reference strand intervals remain an unresolved product class; no intact excised oligo is implied.'))
        g.react(identifier+'_process','INFERRED — illustrative 5-prime resection endpoints',S.DIGESTION,[(source,R.TEMPLATE),(exonuclease,R.CATALYST)],[out,removed],geometry=retention,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Removed strand products and experimental endpoint heterogeneity are unexpanded.'))
        ready.append((out,fragment))
    assembly=CircularCohesiveAssembly(tuple((m.id,f) for m,f in ready),'Explicit neutral target/recipient order; each declared source-style overlap retained once.')
    con=assembly.product('neutral circular reference')
    joint=g.product('joint','INFERRED — paired reference; nicks at both seams',K.MOLECULAR,state=Drawing(tuple(assembly.rows())),chemical_state=con)
    g.react('anneal','INFERRED — pair the declared overlaps',S.HYBRIDIZATION,[(ready[0][0],R.TEMPLATE),(ready[1][0],R.STRUCTURAL)],[joint],geometry=assembly,relations=(M.COHESIVE_ASSEMBLY,))
    polymerase=g.input('polymerase','source Gibson mixture high-fidelity polymerase',K.CATALYST)
    filled=g.product('filled','INFERRED — exact-end toy reference; nicks remain',K.MOLECULAR,state=Drawing(tuple(assembly.rows())),chemical_state=con)
    g.react('fill','INFERRED — polymerase boundary; no gap in this limiting toy',S.EXTENSION,[(joint,R.TEMPLATE),(polymerase,R.CATALYST)],[filled],relations=(M.STRAND_ORDER,),reason='The exact-end reference contains no gap. Variable resection, gap repair and copying errors require separately declared geometry.')
    ligase=g.input('ligase','source Gibson mixture Taq DNA ligase',K.CATALYST)
    sealed=g.product('sealed','INFERRED — intended sealed neutral reference',K.MOLECULAR,state=closed_circle('conditional reference',tuple(con)),chemical_state=con)
    g.react('seal','INFERRED — reference nick sealing',S.NICK_SEALING,[(filled,R.TEMPLATE),(ligase,R.CATALYST)],[sealed],relations=(M.STRAND_ORDER,))
    host=g.input('host','source cloning-host recovery boundary',K.BIOLOGICAL)
    final=g.product('final','INFERRED — conditional recovered reference',K.MOLECULAR,state=closed_circle('neutral reference; no replicon',tuple(con)),chemical_state=con)
    others=g.product('other_products','unresolved backgrounds and alternative assemblies',K.MOLECULAR,state=limitation('Self-joined recipients, residual substrates, alternative assemblies and host-associated material remain unresolved.'))
    g.react('host_recovery','INFERRED — source transformation/recovery boundary',S.MAGIC,[(sealed,R.TEMPLATE),(host,R.TEMPLATE)],[final,others],reason='The neutral example has no replicon or selectable function. Host entry, repair, replication, recovery and success are not modelled; only the intended reference is carried forward.')
    g.mark_terminal(Terminal(others.id,Disposition.UNRESOLVED,'No yield, purity or complete conversion is claimed. Source Results explicitly discusses negative clones.'))
    g.constructed(final,'Conditional constructed-DNA reference; no functional clone or sequencing-library claim.')
    g.declare_reviewed_scope(('parent','guides','cas9','cut','gel_recovery','recover','recipient','exonuclease','target_process','recipient_process','anneal','polymerase','fill','ligase','seal','host','host_recovery','final'),SOURCE+' Exact sequences, prepared-recipient map and endpoints are invented. No genomic target, functional construct or complete substrate distribution is reproduced.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_catch_cloning.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/ncomms9101">Jiang et al.2015</a>. INFERRED — neutral reference DNA, target boundaries and exact resection endpoints. Guides/PAMs and initial cut-end chemistry remain unassigned. This prepared-recipient example models original cloning CATCH; PFGE target-band isolation and nanopore preparation belong to the separate adaptation.'),
        WorkflowSection('Target delimitation → DNA recovery → reference assembly → host boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory reference claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    s=Source(sources_root()/Path(__file__).resolve().parents[1].name/'ncomms9101_PMC4569707.xml.txt','Published CATCH Results and Methods')
    if have(c,s,label='published cloning branch'):
        t=s.path.read_text();c('original explicitly separates cloning from preparative target-band purification','without the need for size selection on PFGE' in t and 'T5 5′–3′ exonuclease, Taq DNA ligase and high-fidelity polymerase' in t,True)
    c.report()
