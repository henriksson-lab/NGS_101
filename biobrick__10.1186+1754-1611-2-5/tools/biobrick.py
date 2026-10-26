"""BioBrick/RFC10 neutral three-fragment reference, source-specific 3A branch."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,Row,revcomp
from plasmid import Enzyme,Cut,find_all
from cohesive_assembly import CohesiveFragment,CircularCohesiveAssembly
from assembly_geometry import StrandRetention,DuplexEndState
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle,unsealed_circular_duplex
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='biobrick.html'
TITLE='BioBrick / RFC10 — illustrative 3A assembly'
SOURCE='Shetty, Endy and Knight2008 doi10.1186/1754-1611-2-5, Results Assembling BioBrick parts/Figure3 and Methods Assembly of BioBrick parts; inspected OpenWetWare noncoding BioBrick prefix/suffix standard, oldid118773. Separately archived RFC10 original/version comparison deferred.'
EXAMPLE='Invented neutral parts and recipient reference, not published constructs or a functional plasmid. The source-specific three-fragment branch uses RFC10 noncoding flanks. Selection markers, replication, experimental recovery and yield are not supplied by the toy.'
PREFIX='GAATTCGCGGCCGCTTCTAGAG'
SUFFIX='TACTAGTAGCGGCCGCTGCAG'
# These explicit recognition/cut declarations select this example, not a generic
# cloning recipe. Negative offsets are measured from the recognition-site end.
E=Enzyme('EcoRI','GAATTC',-5,-1)
X=Enzyme('XbaI','TCTAGA',-5,-1)
SPE=Enzyme('SpeI','ACTAGT',-5,-1)
P=Enzyme('PstI','CTGCAG',-1,-5)

def dna(name,seq):return Construct((Segment(name,seq,inferred=True),),name)

def prepared(identifier,sequence,left,right):
    """Derive physical retained intervals from two locally declared cut sites."""
    parent=dna(identifier+' starting reference',sequence)
    def cuts(enzyme):
        # Explicit local palindromic-site branch, including PstI 3-prime ends.
        if enzyme.site!=revcomp(enzyme.site):raise ValueError('local example requires palindromic recognition sites')
        return [Cut(enzyme.name,p,1,p+len(enzyme.site)+enzyme.top_cut,p+len(enzyme.site)+enzyme.bottom_cut,'') for p in find_all(sequence,enzyme.site,circular=False)]
    lm=cuts(left)
    rm=cuts(right)
    if len(lm)!=1 or len(rm)!=1:raise ValueError('each toy preparation requires exactly the selected two sites')
    lo,hi=lm[0],rm[0]
    if max(lo.top_cut,lo.bottom_cut)>=min(hi.top_cut,hi.bottom_cut):raise ValueError('toy cut order must delimit a nonempty core')
    core=sequence[max(lo.top_cut,lo.bottom_cut):min(hi.top_cut,hi.bottom_cut)]
    fragment=CohesiveFragment(identifier,sequence[min(lo.top_cut,lo.bottom_cut):max(lo.top_cut,lo.bottom_cut)],tuple(dna(identifier+' retained core',core)),sequence[min(hi.top_cut,hi.bottom_cut):max(hi.top_cut,hi.bottom_cut)],'5-prime' if lo.top_cut<lo.bottom_cut else '3-prime','5-prime' if hi.top_cut<hi.bottom_cut else '3-prime',inferred=True)
    geom=StrandRetention(identifier+'_parent',identifier+'_cut',sequence,revcomp(sequence),(lo.top_cut,hi.top_cut),(len(sequence)-hi.bottom_cut,len(sequence)-lo.bottom_cut),fragment,'Explicit invented parent coordinates and local enzyme-cut declarations; not an experimental digestion distribution.')
    return parent,fragment,geom

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    # Recipient reference is opened outside the retained arc; full functional
    # parent maps and the displaced recipient cassette are not fabricated.
    specs=(('prefix',PREFIX+'ACGTCAGTACGATCGT'+SUFFIX,E,SPE),('suffix',PREFIX+'CGATGACCTGATACGA'+SUFFIX,X,P),('recipient',P.site+'GATCGTAGCTCAGTCA'+E.site,P,E))
    ready=[];fragments=[]
    for identifier,sequence,left,right in specs:
        parent,fragment,geom=prepared(identifier,sequence,left,right)
        inp=g.input(identifier+'_parent','invented '+identifier+' parent reference arc',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(parent)),inferred=True,source=EXAMPLE)
        enzymes=g.input(identifier+'_enzymes',left.name+' + '+right.name,K.CATALYST)
        out=g.product(identifier+'_cut','INFERRED — cut '+identifier+' reference',K.MOLECULAR,state=DuplexEndState(fragment,'phosphate','OH','phosphate','OH'),chemical_state=fragment)
        other=g.product(identifier+'_other','other parent intervals and unresolved digest products',K.MOLECULAR,state=Drawing((Row(chunks=[('Unretained parent intervals are not concatenated into an invented excised molecule; complete parent maps and quantitative cleavage fates are unassigned.',None,False)]),)))
        g.react(identifier+'_digest','INFERRED — source enzyme pair on the declared reference arc',S.DIGESTION,[(inp,R.TEMPLATE),(enzymes,R.CATALYST)],[out,other],geometry=geom,relations=(M.STRAND_RETENTION,))
        g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Parent remainders, uncut molecules and alternative digests are unexpanded; no physical purification of these classes is claimed.'))
        if identifier=='recipient':
            phosphatase=g.input('phosphatase','source Antarctic phosphatase',K.CATALYST)
            treated=g.product('recipient_ready','INFERRED — recipient with 5-prime hydroxyls',K.MOLECULAR,state=DuplexEndState(fragment,'OH','OH','OH','OH'),chemical_state=fragment)
            g.react('dephosphorylate','Source destination-vector dephosphorylation',S.DEPHOSPHORYLATION,[(out,R.TEMPLATE),(phosphatase,R.CATALYST)],[treated],reason='Selected published 3A branch; recipient 5-prime hydroxyls cannot supply both phosphates needed to seal both strands at its seams.')
            out=treated
        clean=g.product(identifier+'_purified','INFERRED — selected '+identifier+' class after source purification',K.MOLECULAR,state=out.state,chemical_state=fragment)
        waste=g.product(identifier+'_cleanup_waste','cleanup waste and unassigned residual classes',K.MOLECULAR,state=Drawing((Row(chunks=[('MinElute cleanup partition and residual contaminants are not quantified; no perfect molecular selection is implied.',None,False)]),)))
        g.react(identifier+'_cleanup','Source preparation purification',S.CLEANUP,[(out,R.TEMPLATE)],[clean,waste])
        g.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,'Source cleanup is represented without assigning quantitative retained/discarded fractions.'))
        out=clean
        ready.append(out);fragments.append(fragment)
    geom=CircularCohesiveAssembly(tuple((m.id,f) for m,f in zip(ready,fragments)),'Declared prefix→suffix→recipient order; mixed 5-prime and 3-prime cohesive ends are explicit. Matching ends alone do not establish covalent closure.')
    con=geom.product('intended 3A reference')
    if 'TACTAGAG' not in con.top():raise ValueError('RFC10 noncoding prefix/suffix junction must retain its eight-base scar')
    joint=g.product('joint','INFERRED — annealed three-fragment reference',K.MOLECULAR,state=Drawing(tuple(geom.rows())),chemical_state=con)
    g.react('anneal','INFERRED — pair the selected source-style ends',S.HYBRIDIZATION,[(ready[0],R.TEMPLATE),(ready[1],R.STRUCTURAL),(ready[2],R.STRUCTURAL)],[joint],geometry=geom,relations=(M.COHESIVE_ASSEMBLY,))
    ligase=g.input('ligase','source T4 DNA ligase',K.CATALYST)
    ligated=g.product('ligated','INFERRED — incompletely sealed circular product class',K.MOLECULAR,state=unsealed_circular_duplex(tuple(con),statement='Recipient 5-prime hydroxyls leave at least the corresponding recipient-strand nicks unsealed. Exact ligation fraction and nick coordinates are not asserted.'),chemical_state=con)
    g.react('ligate','INFERRED — source ligation with dephosphorylated recipient',S.NICK_SEALING,[(joint,R.TEMPLATE),(ligase,R.CATALYST)],[ligated],relations=(M.STRAND_ORDER,))
    host=g.input('host','source E. coli transformation/recovery host',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended recovered circular reference',K.MOLECULAR,state=closed_circle('conditional toy circular reference',tuple(con)),chemical_state=con)
    background=g.product('background','alternative assemblies and cell-associated material',K.MOLECULAR,state=Drawing((Row(chunks=[('Cellular repair, selection, parental background and recovery remain unexpanded.',None,False)]),)))
    g.react('recover','INFERRED — host repair and recovery boundary',S.MAGIC,[(ligated,R.TEMPLATE),(host,R.TEMPLATE)],[final,background],reason='Host repair and recovery are unmodelled. The neutral reference lacks a replicon or selection marker; no colony or functional vector is predicted.')
    g.mark_terminal(Terminal(background.id,Disposition.UNRESOLVED,'No purity, clone identity, success rate or complete conversion asserted.'))
    g.constructed(final,'Conditional constructed-DNA reference, with the RFC10 scar preserved; no NGS output.')
    g.declare_reviewed_scope(('prefix_parent','prefix_enzymes','prefix_digest','prefix_cleanup','suffix_parent','suffix_enzymes','suffix_digest','suffix_cleanup','recipient_parent','recipient_enzymes','recipient_digest','recipient_cleanup','phosphatase','dephosphorylate','anneal','ligase','ligate','host','recover','final'),SOURCE+' Exact neutral bodies and recipient reference arc are invented. This is part assembly, not construction of the paper’s functional base vectors.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>','<p class="research-notes"><a href="01_biobrick.html">Research notes</a></p>',info('<a href="https://doi.org/10.1186/1754-1611-2-5">Shetty et al.2008</a>. INFERRED — neutral bodies and recipient reference arc. The RFC10 noncoding scar is retained. Destination dephosphorylation leaves an unsealed product before the host boundary; matching cohesive ends do not validate recovery.'),WorkflowSection('Three source-specific digest pairs → recipient dephosphorylation → ligation → host boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory toy molecular relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    # Reject an extra internal source site instead of silently selecting a cut.
    try:prepared('bad',PREFIX+'GAATTC'+SUFFIX,E,SPE)
    except ValueError:rejected=True
    else:rejected=False
    c('ambiguous parent cleavage is rejected',rejected,True)
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'Biobrick_standard.html.txt','OpenWetWare noncoding BioBrick flanks')
    if have(c,source,label='source standard strings'):
        text=''.join(source.path.read_text().split())
        c('noncoding prefix transcribes the inspected standard',PREFIX in text,True)
        c('suffix transcribes the inspected standard',SUFFIX in text,True)
    c.report()
