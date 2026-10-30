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
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass
from protocol_graph import Terminal,Disposition
TITLE='TRAPS-seq — source reporter interface and unresolved joint-library continuation'
SOURCE='Published SI MOESM1 Fig3caption202–203,SupplementaryMethodsselectedvalidationsections andTable2originalPDFp28;definingpublishedabstract/identityscreened,fullmainMethodsunread'
LIMIT='Private partial source graph. PublishedSIprovidesreportersandvalidationassays;corelibraryorder/pulseassignment/captureattachmentfromdefiningmainremainunread. NoCITE-seq/SEC-seq/kitrecipefilledbyanalogy. ReporterDNA,proteinsecretionandcellRNAareseparatephysicalclasses.'
HANDLE='CCTTGGCACCCGAGAATTCCA'
IFNG1='ACGAGTAGCTCGAGAC'
IFNG2='AGATCAAGAAGCTTCG'
@dataclass(frozen=True)
class DetectionDNA:
    member:int=1
    def __post_init__(self):
        if self.member not in (1,2):raise ValueError('SourceIFNGbarcode1or2required')
    @property
    def core(self):return IFNG1 if self.member==1 else IFNG2
    def scene(self):
        sc=Scene();sc.strand('tag',[Segment('sourcePCRhandle',HANDLE),Segment('sourceIFNGreporteridentity',self.core,feature=feature('traps.ifng_'+str(self.member),'feature_barcode','fixed',note='Source namedbarcode1/2;exacttimepulseassignmentunreviewed,notUMI.')),Segment('sourcepolyA','A'*27)],label='Source5primeBiotin-TEG;'+HANDLE+self.core+'A'*25+'*A*A',mod5='Biotin-TEG')
        sc.note('tag','Two terminalphosphorothioatelinks onsourceDNA;27Aasprinted. Exactantibodyassociationgeometryunreviewed;notdirectprotein→DNAconversion.');sc.physical_reporter=self;return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input('cells','Source secretingcellpopulation;stimulationconditions/sourceexperimentunassigned',K.BIOLOGICAL,primary=True)
    capture=g.input('capture_antibodies','Source cell-surfacecytokinecaptureantibodypreparation;fullattachmentgeometry/coreprocedureunread',K.MOLECULAR)
    trapped=g.product('surface_capture','Sourcecell-surfacecapturedsecretedproteinclasses;freeproteins/bufferpossiblyremain',K.MOLECULAR)
    g.react('capture','Reportedcell-surfacesecretiontrapping;exactcoreoperationsandphysicalseparationsunassembled',S.MAGIC,((cells,R.TEMPLATE),(capture,R.STRUCTURAL)),(trapped,),reason='Definingpublishedabstractsupportsinterface,butfullcoreMethodsunread;SIvalidationcontrolsnotpromotedintotheNGSrecipe.')
    antibody=g.input('detection_antibody','Source antiIFNGdetectionantibody;associationwithtagnotfullyassembled',K.MOLECULAR)
    dna=g.input('detection_dna','PublishedTable2antiIFNGbarcode1DNA',K.OLIGO,state=DetectionDNA().scene())
    sa=g.input('streptavidin','Source streptavidin for antibody coupling;exact reagent/sites unreviewed',K.MOLECULAR)
    absa=g.product('antibody_sa','Source covalent antibody–streptavidin preparation;unbound reactants may remain',K.MOLECULAR)
    g.react('antibody_sa_coupling','PublishedSIFig3: antibody covalently coupled to streptavidin before oligo association',S.MAGIC,((antibody,R.TEMPLATE),(sa,R.STRUCTURAL)),(absa,),reason='Narrow covalent Ab–streptavidin relation stated in published caption;exact coupling reagent/site/workup unreviewed.')
    reporter=g.product('reporter','Source antibody–streptavidin-associated biotinDNA detection preparation;whole complex geometry unassembled',K.MOLECULAR,state=DetectionDNA().scene())
    g.react('reporter_assembly','Source oligo-barcode association after antibody–streptavidin coupling',S.GENERIC,((absa,R.TEMPLATE),(dna,R.STRUCTURAL)),(reporter,))
    probed=g.product('probed','Reporter-probedcapturedsecretion oncells;unboundreporterandbackgroundnotdeclaredremoved',K.MOLECULAR)
    g.react('probe','ReportedDNA-barcodedantibodyprobingofsurfacecapturedproteins',S.MAGIC,((trapped,R.TEMPLATE),(reporter,R.STRUCTURAL)),(probed,),reason='Broadreportedorderonly;sourcecorewash,pulseandtaghandlingoperationsunread,noidealbindingorpurityasserted.')
    library=g.product('joint_output','Reportedjointtranscriptome/reporterlibrarycontinuation;wholemixturesandendgeometryunassembled',K.MOLECULAR)
    g.react('joint_processing','Unreviewedsourcejointsingle-cellbarcoding/copying/librarypreparation',S.MAGIC,((probed,R.TEMPLATE),),(library,),reason='Authoritativemainrecipeunread;SIqPCR/FAM-polyTvalidationnotusedasNGSconstruction. ExactDNA-vsRNAcopyingandlaterPCR/cleanuporderunassembled.')
    g.sequence(library,SequencingPlan(acquisition='unspecified',limitation='Publishedmethodreportssequencingofsecretionsandtranscriptomes;sequencer,cycles,indexes,fullpolymer/read-primergeometryunassembled. SourceTable2libraryprimernamesalonearenotreactionorder.'))
    g.declare_reviewed_scope(('detection_dna','antibody_sa_coupling','reporter_assembly','probe'),SOURCE)
    g.add_source_item(SourceItem('defining_methods_unread','Definingpublishedcorecapture/pulse/libraryMethods','Fullarticle notavailablethroughinitialpublisheraccess;do notfillbyCITE-seq/SEC-seqanalogy.'))
    g.add_source_item(SourceItem('other_sources','Exactlibrary/primer/end geometry andlinkedoriginals','SourceSIvalidationandreporterTable2onlyinspected;nodeclarationofmissinginformation.'))
    return g

def render_page():return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_traps.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'joint_output'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Two source namedIFNGreportersremaindistinct',DetectionDNA(1).core!=DetectionDNA(2).core,True)
    c('ReporterDNA association follows source antibody–streptavidin coupling',g.steps['reporter_assembly'].inputs[0].material,'antibody_sa')
    c('Actual source tag retains5biotinTEG',g.materials['detection_dna'].state.strands['tag'].mod5,'Biotin-TEG')
    c('Source tag identities notUMIs',g.materials['detection_dna'].state.strands['tag'].segs[1].feature.role,'feature_barcode')
    c('No qPCRvalidationstep inserted','q_pcr' in g.steps,False)
    c('Mainunreadpreventsfinishedpublication',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41592-023-01841-y_41592_2023_1841_MOESM1_ESM.pdf.txt'
    if p.exists():
        t=p.read_text();c('Publishedfirstreporterfragmentsverbatim',HANDLE+IFNG1[:11] in t and IFNG1[11:]+'A'*25+'*A*A' in t,True);c('Publishedfirstreportermodificationlabels','5’-Biotin-TEG' in t and '3’-phosphorothioate' in t,True)
    else:c.skip('PublishedSITable2transcription','Externaloriginalabsent')
    c.report()
