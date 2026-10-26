"""Original RF cloning, with duplex megaprimer linear copying kept distinct from PCR."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,Row,revcomp
from pcr import PCRGeometry,PCRPrimer
from assembly_geometry import pcr_product
from protocol_graph import ProtocolGraph,Disposition,Terminal,InputRole as R,MaterialKind as K,StepKind as S,MolecularRelation as M
from protocol_model import Drawing,closed_circle,unsealed_circular_duplex
from protocol_views import WorkflowSection,assemble_page
from page import head,info
PAGE_FILE='rf_cloning.html'
TITLE='RF cloning — illustrative duplex-megaprimer copying'
SOURCE='van den Ent and Lowe2006 doi10.1016/j.jbbm.2005.12.008; published author-hosted original, Methods2.1/2.2 and Results/Figure1, pp68–71. Original RF linear-amplification branch, not EMP or an exponential whole-plasmid PCR.'
EXAMPLE='Invented neutral DNA maps, insert and primers, not the published constructs or an experimental instance. The example lacks a replicon and functional genes. Matching priming windows and intended copying arithmetic do not predict reaction success or recovery.'
LEFT='AGTCGATACGTATCAGCTAGTCGA'
RIGHT='CTAGTCGATGACATCGTACGATGC'
INSERT='ACGTCAGTACGATCGTGACTACGA'
BACKBONE='TCGAGTACCTGACATGTCAGATCA'
DISPLACED='GATCAGTACG'

def dna(name,seq):return Construct((Segment(name,seq,inferred=True),),name)

@dataclass(frozen=True)
class RFReference:
    """Local copying reference only; no amplification kinetics or nick assignment."""
    recipient: Construct
    megaprimer: Construct
    left: str
    insert: str
    right: str
    retained_arc: str
    replaced: str
    def __post_init__(self):
        if self.recipient.top()!=self.left+self.replaced+self.right+self.retained_arc:
            raise ValueError('recipient arc does not match the declared insertion boundaries')
        if self.megaprimer.top()!=self.left+self.insert+self.right:
            raise ValueError('duplex megaprimer does not retain both declared vector priming windows')
        if not all((self.left,self.insert,self.right,self.retained_arc,self.replaced)):
            raise ValueError('RF reference intervals must all be explicit')
        self.scenes # Scene independently checks the two selected 3-prime pairings.
    def reference(self):
        return Construct((Segment('left recipient flank',self.left,inferred=True),Segment('neutral insert',self.insert,inferred=True),Segment('right recipient flank',self.right,inferred=True),Segment('copied recipient arc',self.retained_arc,inferred=True)),'intended RF reference')
    @property
    def scenes(self):
        a=Scene()
        a.strand('recipient bottom',(Segment('right binding bottom',revcomp(self.right),inferred=True),),label='recipient right priming window only')
        a.anneal('top megaprimer',(Segment('left megaprimer tail',self.left,inferred=True),Segment('insert tail',self.insert,inferred=True),Segment('right 3-prime binding',self.right,inferred=True)),to='recipient bottom',pair=('right 3-prime binding','right binding bottom'),label='selected megaprimer strand')
        a.arrow('top megaprimer','copy retained recipient arc; stop boundary unmodelled')
        b=Scene()
        b.strand('recipient top',(Segment('left binding top',self.left,inferred=True),),label='recipient left priming window only')
        b.anneal('bottom megaprimer',(Segment('right megaprimer tail',revcomp(self.right),inferred=True),Segment('insert tail bottom',revcomp(self.insert),inferred=True),Segment('left 3-prime binding',revcomp(self.left),inferred=True)),to='recipient top',pair=('left 3-prime binding','left binding top'),label='other megaprimer strand')
        b.arrow('bottom megaprimer','copy opposite recipient arc; stop boundary unmodelled')
        return a,b
    def rows(self):
        return [row for scene in self.scenes for row in scene.rows()]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,required_output=Disposition.CONSTRUCT,representation='illustrative_example',example_description=EXAMPLE)
    template=dna('neutral insert template',INSERT)
    donor=g.input('donor','invented neutral insert template',K.MOLECULAR,primary=True,state=Scene.duplex(tuple(template)),chemical_state=template,inferred=True,source=EXAMPLE)
    f=dna('forward primer',LEFT+INSERT[:12]);r=dna('reverse primer',revcomp(RIGHT)+revcomp(INSERT[-12:]))
    fm=g.input('forward','invented tailed forward primer',K.OLIGO,state=f,inferred=True,source=EXAMPLE)
    rm=g.input('reverse','invented tailed reverse primer',K.OLIGO,state=r,inferred=True,source=EXAMPLE)
    kod=g.input('kod','source KOD hot-start polymerase for insert PCR',K.CATALYST)
    pg=PCRGeometry.exact(template,PCRPrimer('forward',f.top(),12,inferred=True),PCRPrimer('reverse',r.top(),12,inferred=True),statement='Exact sites in the neutral reference; the toy binding lengths are modelling declarations, not experimental recommendations.')
    amp=pcr_product(pg,f,r,name='duplex megaprimer reference')
    megaprimer=g.product('megaprimer','INFERRED — duplex PCR-product megaprimer',K.MOLECULAR,state=Scene.duplex(tuple(amp)),chemical_state=amp)
    g.react('insert_pcr','INFERRED — tailed insert PCR',S.PCR,[(donor,R.TEMPLATE),(fm,R.STRUCTURAL),(rm,R.STRUCTURAL),(kod,R.CATALYST)],[megaprimer],geometry=pg,relations=(M.PCR_PRODUCT,))
    purified=g.product('purified','INFERRED — purified duplex-megaprimer class',K.MOLECULAR,state=megaprimer.state,chemical_state=amp)
    pcr_waste=g.product('pcr_waste','PCR cleanup waste and unresolved residuals',K.MOLECULAR,state=Drawing((Row(chunks=[('Free primers and cleanup partition are not quantified; no zero-carryover claim.',None,False)]),)))
    g.react('purify','Source PCR-product purification',S.CLEANUP,[(megaprimer,R.TEMPLATE)],[purified,pcr_waste])
    g.mark_terminal(Terminal(pcr_waste.id,Disposition.UNRESOLVED,'Source cleanup exists; exact retained and discarded fractions are unassigned.'))
    recipient_con=dna('recipient reference',LEFT+DISPLACED+RIGHT+BACKBONE)
    recipient=g.input('recipient','invented circular recipient reference; source methylated-parent class',K.MOLECULAR,primary=True,state=closed_circle('parent reference, methylation sites unassigned',tuple(recipient_con)),chemical_state=recipient_con,inferred=True,source=EXAMPLE)
    geom=RFReference(recipient_con,amp,LEFT,INSERT,RIGHT,BACKBONE,DISPLACED)
    active=g.product('active','INFERRED — separate selected megaprimer-strand binding views',K.MOLECULAR,state=Drawing(tuple(geom.rows())),chemical_state=geom)
    g.react('bind','INFERRED — denature and pair selected megaprimer strands',S.HYBRIDIZATION,[(recipient,R.TEMPLATE),(purified,R.STRUCTURAL)],[active],reason='Two panels show separate priming relations, not simultaneous linear placement of an entire circular hybrid. Competing pairings are unassigned.')
    pfu=g.input('pfu','source PfuTurbo for linear amplification',K.CATALYST)
    con=geom.reference()
    copied=g.product('copied','INFERRED — double-nicked circular reference class',K.MOLECULAR,state=unsealed_circular_duplex(tuple(con),statement='Source Figure1 identifies a double-nicked circular duplex. Exact nick coordinates and polymerase stopping are unassigned; no covalent closure is drawn.'),chemical_state=con)
    parent=g.product('parent_carry','methylated parental recipient carryover class',K.MOLECULAR,state=recipient.state,chemical_state=recipient_con)
    residual=g.product('residual','residual megaprimer and alternative copying products',K.MOLECULAR,state=Drawing((Row(chunks=[('Unextended primers, hybrid species and displacement products are unresolved; the replaced parent interval is not an excised molecule.',None,False)]),)))
    g.react('linear_copy','INFERRED — duplex-megaprimer linear amplification',S.EXTENSION,[(active,R.TEMPLATE),(pfu,R.CATALYST)],[copied,parent,residual],geometry=geom,reason='Original RF copies the parental circle using the PCR product as a primer pair. This is not exponential whole-plasmid PCR. Local geometry checks intended base order and two 3-prime binding windows only; enzyme stopping, exact nicks and population kinetics are unmodelled.')
    g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'No yield, quantitative kinetics or complete conversion asserted.'))
    dpni=g.input('dpni','source DpnI',K.CATALYST)
    treated=g.product('treated','INFERRED — copied class after DpnI treatment',K.MOLECULAR,state=copied.state,chemical_state=con)
    background=g.product('background','parental digest, residual parent and hybrid background',K.MOLECULAR,state=Drawing((Row(chunks=[('Methylation sites, DpnI cut products and partition of hemi-methylated hybrids are unassigned.',None,False)]),)))
    g.react('dpni_treatment','Source methylation-dependent parental-template digestion',S.DIGESTION,[(copied,R.TEMPLATE),(parent,R.TEMPLATE),(dpni,R.CATALYST)],[treated,background])
    g.mark_terminal(Terminal(background.id,Disposition.UNRESOLVED,'DpnI treatment is not an asserted physical removal or proof of zero parental background.'))
    host=g.input('host','source E. coli recovery host',K.BIOLOGICAL)
    final=g.product('final','INFERRED — intended recovered circular reference',K.MOLECULAR,state=closed_circle('conditional RF reference',tuple(con)),chemical_state=con)
    other=g.product('other','alternative recovered and cell-associated material',K.MOLECULAR,state=Drawing((Row(chunks=[('Host repair, alternative products and clone recovery are unexpanded.',None,False)]),)))
    g.react('recover','INFERRED — host repair and recovery boundary',S.MAGIC,[(treated,R.TEMPLATE),(host,R.TEMPLATE)],[final,other],reason='The source transforms the nicked product; intracellular repair and recovery are unmodelled. Toy DNA has no replicon or functional phenotype.')
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'No actual clone or functional construct claimed.'))
    g.constructed(final,'Conditional constructed-DNA reference, not an experimental instance or sequencing library.')
    g.declare_reviewed_scope(('donor','forward','reverse','kod','insert_pcr','purify','recipient','bind','pfu','linear_copy','dpni','dpni_treatment','host','recover','final'),SOURCE+' Exact maps and neutral primers are invented. The source 24-base vector windows are retained as a reference choice; copied strand stopping and host chemistry remain unassigned.')
    return g.require_publishable()

def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>','<p class="research-notes"><a href="01_rf_cloning.html">Research notes</a></p>',info('<a href="https://doi.org/10.1016/j.jbbm.2005.12.008">van den Ent and Lowe2006</a>. INFERRED — neutral maps, insert and primers. Original RF uses a duplex megaprimer for linear amplification; exact stopping, nick positions, parental background and host repair remain unassigned.'),WorkflowSection('Insert PCR → duplex-megaprimer linear copying → DpnI → host boundary'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from chemistry_validation import check_graph,Status
    from storage import sources_root
    c=Check();c('no contradictory declared molecular claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    try:RFReference(dna('recipient',LEFT+DISPLACED+RIGHT+BACKBONE),dna('bad megaprimer',LEFT+INSERT+revcomp(RIGHT)),LEFT,INSERT,RIGHT,BACKBONE,DISPLACED)
    except ValueError:rejected=True
    else:rejected=False
    c('wrong megaprimer flank is rejected',rejected,True)
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'rfcloning2006_published.pdf.txt','Published RF Figure1 and Methods2.2')
    if have(c,source,label='source RF linear-copying branch'):
        text=' '.join(source.path.read_text().split())
        c('source distinguishes linear copying and double-nicked output','linear amplification reaction' in text and 'double nicked, circular, double-stranded DNA' in text,True)
    c.report()
