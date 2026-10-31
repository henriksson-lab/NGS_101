from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,extra=()):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,tuple(inputs)+tuple(extra),(y,),changes_molecule=kind is not S.GENERIC,reason='Source-positive operational interface; named supplemental chemistry and complete molecular inventory/geometry remain unreviewed.' if kind is S.MAGIC else '');return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g,protein_present):
    g.declare_reviewed_scope(COMMON_SCOPE+(('protein','equilibrate')if protein_present else ()),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label,kind=S.BINDING):
    y=g.product(key+'_product',label+';unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,kind,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,),reason='Actual main reports RNA-to-DNA copying/amplification; complete primer, enzyme/nucleotide inventory and molecule geometry remain in unread SI.' if kind is S.MAGIC else '');return y

def split(g,key,x,labels):
    out=[g.product(key+'_'+str(i+1),label,K.MOLECULAR) for i,label in enumerate(labels)]
    g.react(key,'Physical source allocation;amounts explicit only where inspected',S.SPLIT,((x,R.TEMPLATE),),tuple(out));return out

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='RNA Bind-n-Seq — associated-manuscript RNA selection interface'
SOURCE='Associated author manuscript PMC4142047 overview54–58,ExperimentalProcedures149–155;actualsupplements01/02/03 gated'
LIMIT='Private incomplete source graph ofselectedCELF1/RNA binding workflow. Associated manuscript version explicit;publishedfinalMethods/fullSIapplicabilityunreviewed. RNA substrate,noncovalentproteinassociation,newRT-cDNA andPCRcopies distinct. No computed affinity,foldingoccupancy,perfectcaptureorpurity. No DNABind-n-Seq/CLIP recipe inheritance.'
ENDPOINT='library_product'
@dataclass(frozen=True)
class SelectedRNA:
    bases:str
    def __post_init__(self):
        if not self.bases or set(self.bases)-set('ACGU'):raise ValueError('ExplicitRNAlettersrequired')
    def conditional_cdna(self):return ''.join({'A':'T','C':'G','G':'C','U':'A'}[b] for b in self.bases[::-1])

def protein_graph(protein_present=True):
    if type(protein_present)is not bool:raise ValueError('Explicit source protein/no-proteincondition required')
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    template=g.input('ivt_template','Source DNAtemplateforrandomRNAinput;exactbody/promoter/endgeometryunreviewed',K.OLIGO,primary=True)
    x=op('ivt',template,'Sourcein-vitrotranscribedrandomRNA input;selected40nt bindingwindow,fullhandles/IVTconditionsdeferred','SourceIVT reagent contributions; complete RNA-input production/workup remains unreviewed',kind=S.MAGIC,reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
    input_control,selected=split(g,'input_allocate',x,('SourceRNAinputcomparisonaliquot;actualpreparation/readoutnotassigned','RNA allocatedtoselectedbindingcondition'))
    g.terminals[input_control.id]=Terminal(input_control.id,Disposition.STORED,'RNAinputcomparison retainedindependently;sourceinputlibrarycontinuationunassembled.')
    buffer=g.input('binding_buffer','Actual source binding-buffer allocation; formula and stockvolume unreviewed',K.MOLECULAR)
    inhibitor=g.input('rnase_inhibitor','Actual source40U Superasin contribution in binding conditions',K.MOLECULAR)
    if protein_present:
        p=g.input('protein','SourceCELF1 purified afterGST-SBPexpression/GSTcleavage;SBPcaptureassociationnoncovalent',K.MOLECULAR)
        p=op('equilibrate',p,'Sourceproteinconditionequilibrated30minRT;concentrationunassigned',extra=((buffer,R.STRUCTURAL),))
        x=join(g,'bind',(selected,p,inhibitor),'SourceRNA1µM+40USuperasin,1hRT;sourceCELF1route,notRBFOX37°Ccondition')
    else:x=op('bind',selected,'Sourcezeroproteincontrolcondition;RNA1µM/1hRT plusbuffer/inhibitor,notfictionalproteininput',extra=((inhibitor,R.STRUCTURAL),(buffer,R.STRUCTURAL)))
    beads=g.input('streptavidin','Source1mgwashedSA magneticbeads',K.MOLECULAR)
    x=join(g,'capture',(x,beads),'Source1hSBPprotein/RNAcapture or zero-proteinbackgroundcondition')
    x=clean('wash',x,'Sourcecollectbeads/removeunboundRNA andone1mLwash;retainbeadpreparation')
    hc=g.input('heat_condition','Source70C10min elution conditions',K.CONDITION)
    x=op('heat',x,'Source70°C10min100µL10mMTrispH7/1mMEDTA/1%SDSelution;beads/denaturedprotein/reagentsremainuntilcollection','Actual100ulTris/EDTA/SDS elution contribution',kind=S.THERMAL_TREATMENT,reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL,extra=((hc,R.CONDITION),))
    y=g.product('eluate','Sourcecollectedelutedmaterial containingRNA andpossibleprotein/remnants;notyetextracted',K.MOLECULAR);rest=g.product('bead_remainder','Residualbeads/liquid/unreleasedmaterial;sourcecollectionseparatesfromeluate',K.MOLECULAR)
    g.react('collect','Source physically collectelutedmaterial',S.CLEANUP,((x,R.TEMPLATE),),(y,rest));g.terminals[rest.id]=Terminal(rest.id,Disposition.UNRESOLVED,'Residualfractionnotautomaticallydestroyedorsequenced.')
    x=clean('rna_extract',y,'SourceboundRNA extraction;exactreagents/partitiongeometryunreviewed')
    rtprimer=g.input('rt_primer','SourceRTprimerexactbodyunreviewed',K.OLIGO)
    x=join(g,'rt',(x,rtprimer),'SourceRTcreatesnewDNA complementary toselectedRNA;RNA/protein/remnantsfatenotimplicitlyerased',kind=S.MAGIC)
    f=g.input('pcr_forward','SourcePCRforwardprimerunreviewed',K.OLIGO);r=g.input('pcr_reverse','SourcePCRreverseprimerunreviewed',K.OLIGO)
    x=join(g,'pcr',(x,f,r),'SourcePCRamplifiescDNA;newDNAcopies distinct fromnativeRNA/SBPprotein',kind=S.MAGIC)
    x=op('library',x,'Unreviewedsourcepost-PCRworkup/library/acquisitioncontinuation',kind=S.MAGIC)
    end(g,x,'SourceselectedRNAcDNAsequencedforbindingmotifs;exactplatform/runlength/indices/primerlandingunreviewed. Relativeaffinity/foldingmotifanalysisiscomputational,notphysicalC-to-Tconversion.')
    return incomplete(g,protein_present)

# Literal reviewed targets from associated main153-155; neither opaque copying
# nor known high-level preparation covers the original unassembled source item.
COMMON_SCOPE=tuple('ivt_template ivt_reagent ivt input_allocate binding_buffer rnase_inhibitor bind streptavidin capture wash_selection wash heat_condition heat_reagent heat collect rna_extract_selection rna_extract rt_primer rt pcr_forward pcr_reverse pcr library library_product'.split())

def protocol_graph():return protein_graph()
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();zero=protein_graph(False);r=SelectedRNA('AUGC')
    c('RNA/protein association is binding',g.steps['bind'].kind,S.BINDING)
    c('No-protein control receives actual inhibitor','rnase_inhibitor' in zero.materials,True)
    c('Thermal elution has physical buffer contribution',g.materials['heat_reagent'].kind,K.MOLECULAR)
    c('Unknown source document is not reaction reagent','library_reagent' in g.materials,False)
    c('ConditionalDNAcopyreportscomplementarysequence',r.conditional_cdna(),'GCAT')
    c('CopycalculationretainsoriginalRNA',r.bases,'AUGC')
    try:SelectedRNA('ATGC');bad=False
    except ValueError:bad=True
    c('RNAconstructorrejectsDNAthymidine',bad,True)
    c('Zero-proteincontroldoesnotinventprotein','protein' in zero.materials,False)
    c('Protein capture follows RNA binding',g.steps['capture'].inputs[0].material,'bind_product')
    c('Extractionusesactualcollectedeluate',g.steps['rna_extract'].inputs[0].material,'eluate')
    c('Sourceinputcomparison retained',g.terminals['input_allocate_1'].disposition,Disposition.STORED)
    c('Source-limited geometry prevents promotion',g.publishable,False);c.report()
