"""IVA subcloning: two copying relations in one source PCR mixture."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Row,Scene,Segment,revcomp
from pcr import PCRPrimer,PCRGeometry
from assembly_geometry import pcr_product,CircularOverlapReference
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='iva_cloning.html'
TITLE='IVA cloning — illustrative single-mixture subcloning'
SOURCE='Garcia-Nafria, Watson and Greger2016 doi10.1038/srep27459; Method overview, Sub-cloning/Figure3e, Methods PCR reactions and time courses. Inverse-PCR mutation branches and alternative vector entry routes excluded.'
EXAMPLE='Invented short neutral template maps, primers and overlaps, not the published constructs or an experimental instance. Parent methylation is a source-qualified class precondition; exact methylated sites are unassigned. The toy has no replicon or claimed biological function.'
LEFT='AGTCGATACGTA'
RIGHT='CTAGTCGATGAC'
INSERT='ACGTCAGTACGATCGTGACTACGA'
RECIPIENT='TCGAGTACCTGACATGTCAGATCA'

def dna(name,seq):return Construct((Segment(name,seq,inferred=True),),name)

@dataclass(frozen=True)
class ParentReference:
    con: Construct
    methylation_scope: str
    def __post_init__(self):
        if not self.methylation_scope.strip():raise ValueError('parent methylation scope must be explicit')
    def rows(self):
        return [Row(chunks=[(self.methylation_scope,None,False)]),*closed_circle('toy parental reference',tuple(self.con)).rows()]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    pol=g.input('polymerase','source Phusion high-fidelity DNA polymerase',K.CATALYST)
    amplicons=[];parents=[]
    for identifier,sequence,ftail,rtail in (('insert',INSERT,LEFT,revcomp(RIGHT)),('recipient',RIGHT+RECIPIENT+LEFT,'','')):
        con=dna(identifier+' parent map',sequence)
        parent_state=ParentReference(con,'Source methylated-parent class; exact methylation pattern and complete digestion are unassigned in this reference example.')
        parent=g.input(identifier+'_parent','invented parental reference map',K.MOLECULAR,primary=True,state=parent_state,chemical_state=con,inferred=True,source=EXAMPLE)
        f=dna(identifier+' forward',ftail+sequence[:12]);r=dna(identifier+' reverse',rtail+revcomp(sequence[-12:]))
        fm=g.input(identifier+'_f','invented '+identifier+' primer',K.OLIGO,state=f,inferred=True,source=EXAMPLE)
        rm=g.input(identifier+'_r','invented '+identifier+' primer',K.OLIGO,state=r,inferred=True,source=EXAMPLE)
        geom=PCRGeometry.exact(con,PCRPrimer(f.name,f.top(),12,inferred=True),PCRPrimer(r.name,r.top(),12,inferred=True),statement='Reference map opened on the selected copying arc; source PCR mixture contains both template/primer pairs. Circular template topology, cross-priming and yield are not predicted.',amplicon_id=identifier+'_amp')
        product=pcr_product(geom,f,r,name=identifier+' amplicon reference')
        amp=g.product(identifier+'_amp','INFERRED — '+identifier+' amplicon reference',K.MOLECULAR,state=Scene.duplex(tuple(product)),chemical_state=product)
        carry=g.product(identifier+'_carry','parental-template carryover class',K.MOLECULAR,state=parent_state,chemical_state=con)
        g.react(identifier+'_pcr','INFERRED — '+identifier+' copying relation in the shared PCR mixture',S.PCR,[(parent,R.TEMPLATE),(fm,R.STRUCTURAL),(rm,R.STRUCTURAL),(pol,R.CATALYST)],[amp,carry],geometry=geom,relations=(M.PCR_PRODUCT,))
        amplicons.append(amp);parents.append(carry)
    dpni=g.input('dpni','source DpnI',K.CATALYST)
    treated=[]
    for identifier,amp in zip(('insert','recipient'),amplicons):
        treated.append(g.product(identifier+'_treated','INFERRED — amplicon class after source DpnI treatment',K.MOLECULAR,state=amp.state,chemical_state=amp.chemical_state))
    background=g.product('parent_background','digested and residual parental-template classes',K.MOLECULAR,state=Drawing((Row(chunks=[('Parent fragments, undigested background and quantitative partition are unassigned.',None,False)]),)))
    g.react('dpni_treatment','DpnI treatment of the shared reaction mixture',S.DIGESTION,[(m,R.TEMPLATE) for m in amplicons+parents]+[(dpni,R.CATALYST)],[*treated,background])
    g.mark_terminal(Terminal(background.id,Disposition.UNRESOLVED,'No physical cleanup or zero parental background inferred; methylation-specific cleavage is not base-resolved here.'))
    refs=tuple((m.id,m.chemical_state,len(LEFT),len(RIGHT)) for m in treated)
    # Recipient has swapped end windows; these overlap lengths happen to be equal,
    # while sequence identity is checked independently in the reference geometry.
    geom=CircularOverlapReference(refs,'Invented insert-only tails match the retained recipient ends; source intracellular recombination is unmodelled.')
    con=geom.reference('intended subcloning reference')
    host=g.input('host','source bacterial in-vivo assembly/recovery host',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended assembled circular reference',K.MOLECULAR,state=closed_circle('intended toy reference',tuple(con)),chemical_state=con)
    other=g.product('other','alternative assembly and cell-associated material',K.MOLECULAR,state=Drawing((Row(chunks=[('Alternative molecular fates and intracellular material remain unexpanded.',None,False)]),)))
    g.react('host_assembly','INFERRED — bacterial assembly and recovery boundary',S.MAGIC,[(m,R.TEMPLATE) for m in treated]+[(host,R.TEMPLATE)],[final,other],geometry=geom,reason='Source in-vivo assembly is not an in-vitro ligation. Resection, strand exchange, repair and clone selection remain unmodelled; matching overlaps do not establish recovery. The toy has no replicon.')
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'No quantitative yield, complete conversion or defined off-product frequency claimed.'))
    g.constructed(final,'Conditional reference sequence only; no verified clone, actual construct function or NGS readiness.')
    g.declare_reviewed_scope(('polymerase','insert_parent','insert_f','insert_r','insert_pcr','recipient_parent','recipient_f','recipient_r','recipient_pcr','dpni','dpni_treatment','host','host_assembly','final'),SOURCE+' Exact maps and primer pairs are invented. The two PCR panels are simultaneous copying relations within one source reaction, not sequential/separate PCRs.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_iva_cloning.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/srep27459">Garcia-Nafria et al.2016</a>. INFERRED — invented maps, primers and overlap values. The two PCR panels describe simultaneous reference-copying relations in one tube; they are not separate PCRs. Host chemistry, methylation-specific cleavage, terminal heterogeneity and background frequencies are unassigned.'),
        WorkflowSection('Shared PCR → DpnI mixture treatment → in-vivo assembly boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory example reference relations',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'srep27459_PMC4893743.xml.txt','Published Sub-cloning and PCR methods')
    if have(c,source,label='source shared-reaction route'):
        text=' '.join(source.path.read_text().split())
        c('source explicitly puts both subcloning templates in one PCR','amplifying both insert and vector from two plasmids' in text and 'in a single PCR' in text,True)
        c('source PCR polymerase matches selected branch','Phusion® High-Fidelity DNA polymerase' in text,True)
    c.report()
