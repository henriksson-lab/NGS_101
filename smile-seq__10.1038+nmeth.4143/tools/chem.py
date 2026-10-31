from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,),changes_molecule=False);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(SI_CAPTURE_SCOPE,SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label):
    y=g.product(key+'_product',label+'; unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,S.MAGIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,),changes_molecule=False,reason='Source-reported physical surface/complex preparation; exact immobilization and association chemistry unassembled. No DNA-copy or ligation claim.');return y

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='SMiLE-seq — published microfluidic TF-binding front end'
SOURCE='Published supplement MOESM206 Fig1 setup and Fig5a four-step heterodimer capture caption; main unavailable'
LIMIT='Private incomplete graph of the published heterodimer front end. Exact surface chemistry, DNA polymers/Cy5 attachment position, release and library/NGS operations remain unreviewed. No perfect dimer formation, affinity, membrane trapping or purity is predicted. This TF-binding method is not the unrelated PacBio WGS homonym Smile-seq.'
ENDPOINT='library'
SI_CAPTURE_SCOPE=('antibody','chip','surface','dna_bait','tf1','tf2','load','trap_reagent','trap','observe','wash_selection','wash')
@dataclass(frozen=True)
class ReporterInventory:
    body:str
    reporter:str
    attachment:str='unreviewed'
    def __post_init__(self):
        if (self.body,self.reporter) not in {('dsDNA','Cy5'),('TF1','eGFP'),('TF2','mCherry'),('anti-eGFP antibody','biotin')}:raise ValueError('Selected source reporter/body pair required')
    @property
    def dna_bases(self):return None

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    antibody=g.input('antibody','Source biotinylated anti-eGFP antibody; biotin linkage/surface geometry unreviewed',K.MOLECULAR,chemical_state=ReporterInventory('anti-eGFP antibody','biotin'))
    chip=g.input('chip','Source PDMS chip on plasma-activated glass, with flexible membrane/button and control/flow channels',K.MOLECULAR)
    x=join(g,'surface',(chip,antibody),'Fig5a Step1: source antibody immobilization under button; exact chemistry/workup unassembled')
    dna=g.input('dna_bait','Source Cy5-labelled duplex DNA binding bait; whole sequence and dye position unreviewed, not a UMI/barcode',K.OLIGO,primary=True,chemical_state=ReporterInventory('dsDNA','Cy5'))
    tf1=g.input('tf1','Source eGFP-tagged dimerizing TF1; full protein/association not instantiated',K.MOLECULAR,chemical_state=ReporterInventory('TF1','eGFP'))
    tf2=g.input('tf2','Source mCherry-tagged TF2 partner; full protein/association not instantiated',K.MOLECULAR,chemical_state=ReporterInventory('TF2','mCherry'))
    x=join(g,'load',(x,dna,tf1,tf2),'Fig5a Step2: introduce tagged partners and labelled DNA; bound/unbound species coexist')
    x=op('trap',x,'Fig5a Step3: source complexes trapped under flexible PDMS membrane; selected surviving association class only','Source control microvalves/button; no occupancy or guaranteed co-binding prediction')
    x=op('observe',x,'Source fluorescent on-chip dimer readout; observation does not ligate DNA or alter reporter identity')
    x=clean('wash',x,'Fig5a Step4: wash outside-button unbound molecules and molecular complexes; retain source membrane-protected preparation')
    y=g.product('library','Reported DNA sequencing continuation; physical release/copying/fraction and oligo geometry unassembled',K.MOLECULAR)
    g.react('library_boundary','Unreviewed defining-main release/library/acquisition continuation',S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Published SI establishes real capture operations; full main library procedure unread, no HT-SELEX/Bind-n-Seq/modern kit recipe supplied by analogy.')
    end(g,y,'Published method reports TF-associated DNA sequencing; exact platform, read cycles, indexes, primers, recovery and library handling unreviewed.')
    return incomplete(g)
def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Cy5 reporter is not converted into DNA bases',ReporterInventory('dsDNA','Cy5').dna_bases,None)
    try:ReporterInventory('DNA','eGFP');bad=False
    except ValueError:bad=True
    c('Reject source-incompatible DNA/protein reporter pair',bad,True)
    c('Source surface is prepared before ligand/partner loading',g.steps['load'].inputs[0].material,'surface_product')
    c('Physical washing follows membrane trapping',g.steps['wash'].inputs[0].material,'observe_product')
    c('Library continuation consumes retained fraction',g.steps['library_boundary'].inputs[0].material,'wash_kept')
    c('Unknown dye site is not assigned 5prime',g.materials['dna_bait'].chemical_state.attachment,'unreviewed')
    c('Fluorescent observation does not produce a copied DNA molecule',g.steps['observe'].changes_molecule,False)
    c('Source-limited graph remains private',g.publishable,False);c.report()
