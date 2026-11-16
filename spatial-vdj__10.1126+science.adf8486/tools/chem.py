from pathlib import Path
import html
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'10x-visium-spatial-gene-expression'/'tools'))
import visium_rev_a
from chemdraw import (Scene,Segment,Construct,MolecularState,Workflow,workflow_panel,
                      feature,complement_segments,revcomp)
from protocol_model import role
from protocol_handoff import import_material
from pcr import PCRGeometry
import oligos_spatial_vdj as oligos
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition,AnnotationSlot,SourceClaim,SourceDiscrepancy
from page import head, info, caveat, table
from protocol_views import assemble_page, WorkflowSection, PlanLimitationSection, HandoffSection
from chemdraw import panel
PAGE_FILE = 'spatial-vdj.html'

# Curator-declared source inventory. This is intentionally independent of graph
# enumeration: adding a new step/input requires another located source review.
# Final published Methods: capture1625–1642, SR-TRB1645–1668.
MAIN_COMMON_TARGETS = (
    'source','select10ul','preamp_allocate','preamp_cdna_primers',
    'preamp_kapastock','preamp_components','preamp_programme','preamp',
    'preamp_pool','preamp_clean','preamp_clean_beads','preamp_clean_condition',
    'concentrate','capture_concentration','capture','capture_ratio','xgen_capture',
    'postcapture_cdna_primers','postcapture_kapastock','postcapture_components',
    'postcapture_programme','postcapture','postcapture_clean',
    'postcapture_clean_beads','postcapture_clean_condition','outer_calibration',
    'outer_allocate','pool_replicates','index_kapastock','index_components',
    'index_primers','index_programme','index','index_clean',
    'index_clean_beads','index_clean_condition','final_qc','load_select',
    'other_samples','sequencing_pool','pool_loading','load','library',
)
# Four actual independently processed replicates, as declared by the source.
SOURCE_REPLICATES = (1, 2, 3, 4)
REPLICATE_MAIN_TARGETS = (
    'outer_kapastock','outer_components','outer_programme','outer',
    'outer_clean','outer_clean_beads','outer_clean_condition','outer_qc',
    'inner_calibration','select5ng','inner_kapastock','inner_components',
    'inner_programme','inner','inner_clean','inner_clean_beads',
    'inner_clean_condition','inner_qc','select10ng',
)
# Data S1 supplies these complete reagent patterns, not full template occupancy.
DATA_S1_COMMON_TARGETS = ('ig_probes','tr_probes','capture_blockers')
DATA_S1_REPLICATE_TARGETS = ('outer_v','outer_partread1','inner_v','inner_partread1')

class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL,source_factory=None):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance')
        self.x=source_factory(self.g)if source_factory else self.g.input('source',start,kind,primary=True,chemical_state={})
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

AUTHORITY='Final published Science full article, DOI10.1126/science.adf8486, authorized browser HTML science_adf8486_final_browser.html.txt: Results1580/Methods1623–1669; supplementary index1767–1782. Final published LR processing1685 additionally verifies16/12 analysis metadata; this does not establish full SR primer geometry. Published SM PDF title/Methods1–72 and source DataS1 spatial worksheets inspected after user supply; full kit/template geometry unresolved. Earlier accepted-only review is historical, not whole-version equivalence.'
PREAMP_PROGRAMME='Source literal: Repeat steps 2-5 six times for a total of five cycles. Executed count unassigned; no silent correction.'
def unknown_oligo(label,identifier=None):
    s=Scene();s.strand('source',[Segment(label,'[wholeorderingunassembled]',placeholder=True,feature=identifier)]);return s

def pcr_supplies(b,key,label,primers,kapa=True):
    g=b.g;stock=g.input(key+'_kapastock',('Actual KAPA HiFi HotStart ReadyMix2X KK2602 '+label)if kapa else('Actual indexPCR aggregate enzyme/NT/buffer supplies; identity/packaging/stockcounts unassigned '+label),K.MOLECULAR,chemical_state={'physical_stock':'One actual ReadyMix allocation per reaction; components copresent, no fresh separate additions'if kapa else'Aggregate actual indexPCR supplies; premix status/stockcounts/identities not assigned','volume':'25ul in source common pre/postcapture PCR; outer/inner stock volumes unassigned'if kapa else'IndexPCR supply volumes unassigned','packaging':'Named ReadyMix only, not assumed standalone buffer/enzyme additions'if kapa else'Not claimed to be a KAPAReadyMix or another namedkit'})
    if key=='preamp':stock.chemical_state.update({'physical_stock':'Aggregate actual ReadyMix allocations across a variable collection of separate PCR reactions; each actual reaction has one ReadyMix allocation, not one shared 25ul stock for all templates','volume':'25ulReadyMix PER reaction; multiplicity and aggregate totalamountunassigned','source_array_boundary':'Collection-mode source interface hides separate PCR wells; no duplication of one reaction, no physical pooling beforePCR'})
    enzyme=g.product(key+'_enzyme','KAPA DNApolymerase component of actual same ReadyMix'if kapa else'ActualindexPCRDNApolymeraseclass, identityunassigned',K.CATALYST);nt=g.product(key+'_nt','Nucleotide components of actual same ReadyMix; concentrations unassembled'if kapa else'ActualindexPCRnucleotideclass; concentrationsunassigned',K.MOLECULAR);buf=g.product(key+'_buffer','Buffer components of actual same ReadyMix; formulation unassembled'if kapa else'ActualindexPCRbufferconditions; formulationunassigned',K.CONDITION)
    g.react(key+'_components','Copresent components of actual named ReadyMix allocation; no physical dissociation'if kapa else'ExposeactualindexPCRsupplyclasses; packagingandnumberofstocksunassigned',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nt,buf),changes_molecule=False)
    return ((enzyme,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION))+tuple(primers)
def primer(b,key,label,identifier=None):
    if key in ('preamp_cdna_primers','postcapture_cdna_primers'):
        return b.stock(key,'Visium Rev A cDNA Primers PN-2000089: Partial Read1 forward and Partial TSO reverse',K.OLIGO,R.STRUCTURAL,
            state=visium_rev_a.cdna_primer_scene(),source=visium_rev_a.SOURCE,
            chemical_state={'source_scope':'Same cDNA primer preparation in the exact Visium revision cited by the final Spatial VDJ Methods; whole sequence patterns, not an independent stock-dose declaration','end_modifications':'unassigned'})
    kind = ('partread1' if key.endswith('_partread1') else
            'outer' if key.endswith('_outer_v') else
            'inner' if key.endswith('_inner_v') else
            {'ig_probes':'ig','tr_probes':'tr','capture_blockers':'blockers'}.get(key))
    if kind:
        return b.stock(key,'Published Data S1 ' + kind + ' oligo inventory; end modifications/occupancy unassigned',
                       K.OLIGO,R.STRUCTURAL,state=oligos.pool_scene(kind),
                       chemical_state={'source':'Published Data S1 worksheet whole sequence fields',
                                       'physical_scope':'Actual pool allocation; source inventory rows, including repeated sequences; not a connected molecule or claim of distinct species/molarity',
                                       'end_chemistry':'Not assigned from sequence-only cells'})
    return b.stock(key,label,K.OLIGO,R.STRUCTURAL,state=unknown_oligo(label,identifier))
def amp(b,key,label,primers,programme,kapa=True):
    supplies=pcr_supplies(b,key,'for '+label,primers,kapa=kapa)+(b.cond(key+'_programme',programme),)
    if key.endswith(('_outer','_inner')):
        b.op(key,label,S.PCR,supplies,geometry=PCRGeometry.unavailable(
            unknown_oligo('Actual captured transcript template sequence unassembled'), reason='Published Data S1 gives both primer pools; actual captured transcripts, occupied sites and complete amplicons remain unassembled.'))
    else:
        reason=('The Rev A cDNA primer pair is now declared, but actual source reaction collections, first-cycle occupancy and products remain unassembled; knowing the primer orders does not validate the PCR.'if key in ('preamp','postcapture')else'Published indexing order; selected whole index primers and occupied products remain unassembled. Index supplies are not assumed to be the earlier ReadyMix.')
        b.op(key,label,S.MAGIC,supplies,reason=reason)

def cleanup(b,key,label):b.clean(key,label,extra=(b.stock(key+'_beads','Actualsource '+label+' bead allocation',K.SELECTION,R.SELECTION),))
def qc(b,key,label):
    g=b.g;y=g.product(key+'_analytical',label,K.QC);r=g.product(key+'_stock','Actual preparation remaining after analytical allocation; amountunassigned',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(y.id,Disposition.QC,'Actualanalyticalfraction excludedfromproduction'));b.x=r

def visium_input(graph):
    return import_material(graph,'source',visium_rev_a.full_length_export(),
        expected_protocol=visium_rev_a.PROTOCOL,expected_version='CG000239 Rev A',
        source='Final published Science adf8486 Methods, Preparation of Spatial GEX and target enrichment: CG000239 Rev A;10ul full-length Visium cDNA before fragmentation.',
        label='Full-length Visium Rev A cDNA from upstream Step4.3 cleanup, before fragmentation',primary=True,
        requirements={'stage':'amplified_full_length_cdna','fragmented':False,'spatial_barcode_nt':16,'umi_nt':12,'sample_index_present':False},
        required_sites={'PartRead1':oligos.PARTREAD1,'cDNA reverse':visium_rev_a.REVERSE})


def protocol_graph():
    b=B('Spatial VDJ published SR-TRB four-replicate partial preparation',AUTHORITY,'Full-length Visium Rev A cDNA',K.MOLECULAR,source_factory=visium_input);g=b.g
    b.split('select10ul','Actual10ulfromVisiumfull-lengthcDNAstock','Actualselected10ulcDNA','UntransferredVisiumcDNAstock')
    collection=g.product('preamp_collection','Actual collection of separate10ngtemplate PCR subsamples; number depends on measured concentration, not a physical pool',K.MOLECULAR);g.react('preamp_allocate','Splitselected10ulinto10ngtemplates, multiplicityunassigned',S.GENERIC,((b.x,R.TEMPLATE),),(collection,),changes_molecule=False);b.x=collection
    amp(b,'preamp','Commonpreamplification of separate10ngsubsamples, notpooledbeforePCR',(primer(b,'preamp_cdna_primers','Aggregateactual10xcDNAprimerallocations7.5ulPERPCR acrossvariablearray; number/totalamountunassembled; wholeRevAcDNAprimerordersdeclared'),),PREAMP_PROGRAMME+' Source25ulReadyMix+7.5ulprimers+17.5ulsample/water;98C3min,98C15s/63C30s/72C2min;72C1minfinal.')
    b.op('preamp_pool','Source physical pooling of completedpreamplification subsamples',S.GENERIC,product='ActualpooledpreamplifiedspatialcDNA; parallelreactionmultiplicity/wholeproductsunassembled')
    cleanup(b,'preamp_clean','0.6xSPRIselect; elute40ulQiagenEB')
    b.op('concentrate','SourceIDTxGenv4optionalAppendixAAMPureconcentrationthenPlateProtocolstep13',S.MAGIC,(b.stock('capture_concentration','ActualAMPureconcentrationpreparation; delegatedamounts/completeconditionsunassembled',K.SELECTION,R.SELECTION),))
    b.op('capture','Source16hovernightIGandTRcohybridization/capture; selectedIG:TR1:3',S.MAGIC,(primer(b,'ig_probes','ActualIGconstantregionDiscoverypool; DataS1 whole pool now declared; end modifications unassembled'),primer(b,'tr_probes','ActualTRA/TRBconstantregionDiscoverypool; TRG/TRDnotincluded, DataS1 whole pool now declared; end modifications unassembled'),primer(b,'capture_blockers','Actualcustomuniversal-adapterblockingoligos; DataS1 whole pool now declared; end modifications unassembled'),b.stock('xgen_capture','ActualIDTxGenHybridizationandWashKit1080584 v4capture/washpreparation; selection/packagingcomponentsunassembled',K.SELECTION,R.SELECTION),b.cond('capture_ratio','TRbaitsequencestriplicatedinorder; equalvolumesofIG/TRproducedpoolsachieve1:3, NOTthreefoldTRpoolvolume.')),
         reason='Published-source hybridization capture established; baitmodifications/selectionligands/wholeboundRNA-derivedDNAgeometry and linkedIDTv4recipeunassembled. No assumed biotin/streptavidin or fresh vendor-default components.')
    amp(b,'postcapture','Selectedsecondexperiment:one14cyclepostcapturePCR',(primer(b,'postcapture_cdna_primers','Actual7.5ulVisiumRevAPN2000089cDNAprimerpreparationperPCR; wholeordersdeclared'),),'Source14cyclesbranch; firstexperiment7+5cyclebranchnotmerged; source25ulKAPAReadyMix/17.5ulsampleandwater, programmeotherwiseasreported.')
    cleanup(b,'postcapture_clean','SourceAMPurepostcapturePCRpurification; exactratio/volumesunassembled')
    qc(b,'outer_calibration','ActualqPCRcyclecalibrationanalyticalallocation, amountandcompleteprogrammeunassigned')
    starts=tuple(g.product('outer_template_'+str(i),'Actual2.5–5ngcapturedcDNAallocation forTRBreplicate'+str(i+1),K.MOLECULAR)for i in range(4));rest=g.product('outer_template_rest','UntransferredcapturedcDNAstock',K.MOLECULAR)
    g.react('outer_allocate','ActualfourSRTRBreplicatesfromcapturedstock',S.SPLIT,((b.x,R.TEMPLATE),),starts+(rest,));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredcapturedstock; yieldunassigned'));indexed=[]
    for i,start in enumerate(starts):
        k='rep'+str(i+1);b.x=start
        amp(b,k+'_outer','TRBouterPCR12cyclesreplicate'+str(i+1),(primer(b,k+'_outer_v','ActualOuterTRBVprimerpool5primeofCDR3; DataS1 whole pool now declared; occupied transcript sites unassembled'),primer(b,k+'_outer_partread1','Published whole PartRead1; occupied transcript site unassembled')),'Sourceprintedprimerdilution100uMstocks→40xdilution→2.5uM; finalPCRconcentrationconventionunassigned;98C5min;12x98C20s/65C30s/72C90s;72C7min; TRA14cyclesoutsidechosenTRBscope.')
        cleanup(b,k+'_outer_clean','0.6xAMPure/two80percentEtOH/EB15min37C')
        qc(b,k+'_outer_qc','SourceQubit/Bioanalyzeranalyticalallocation; amountunassigned')
        qc(b,k+'_inner_calibration','SourceinnerPCRqPCRcalibrationanalyticalallocation; amountunassigned')
        b.split(k+'_select5ng','Actual5ngouterproductforinnerPCR','Actual5nginnertemplate','UntransferredouterPCRstock')
        amp(b,k+'_inner','SourceinnerTRBPCReightcycles',(primer(b,k+'_inner_v','ActualInnerTRBVprimerpooladjacenttoCDR3withTruSeq-compatiblehandle; DataS1 whole pool now declared; occupied transcript sites unassembled'),primer(b,k+'_inner_partread1','Published same whole PartRead1; occupied transcript site unassembled')),'Sourceconcentration/reagentsasouter;98C5min;8x98C20s/72C30s/72C90s;72C7min.')
        cleanup(b,k+'_inner_clean','SourceAMPureandEtOHwashasdescribedforouter; exactrecoveryunassigned')
        qc(b,k+'_inner_qc','SourceQubit/Bioanalyzeranalyticalallocation; amountunassigned')
        b.split(k+'_select10ng','Actual10ngfromreplicateforindexingpool','Actual10ngselectedinnerPCRproduct','UntransferredinnerPCRstock');indexed.append(b.x)
    pool=g.product('index_pool','Actualfour10ngreplicateallocationspooledforonesample/sameindex',K.MOLECULAR);g.react('pool_replicates','Sourcepoolbeforeindexing, notfourindependentlyindexedlibraries',S.POOL,tuple((x,R.TEMPLATE)for x in indexed),(pool,));b.x=pool
    amp(b,'index','SourceTruSeqindexPCR8cycles',(primer(b,'index_primers','ActualTruSeqindexingprimerpreparation; selected whole index orders/paired index geometry unassembled',feature('spatialvdj.sample','sample_index','unknown',note='Sourceassignssameindextoallfourreplicates; wholeindexposition/bodyunassembled, notinferredfromTruSeqname')),),'Source98C3min;8x98C20s/60C30s/72C90s;72C5min; completehistoricalindexsupplyformulationunassembled.',kapa=False)
    cleanup(b,'index_clean','SourceAMPureandEtOHwashafterindexPCR')
    qc(b,'final_qc','SourceQubit/Bioanalyzerfinalanalyticalallocation; amountunassigned')
    b.split('load_select','Actualsamplesselectedforsequencingpool; amountsunassigned','SelectedindexedSRTRBsamplelibraryallocation','UntransferredindexedSRTRBlibrary')
    b.op('sequencing_pool','Sourcepoolsindexed samples forsequencing; otheractualsamples/molaritiesunassembled',S.MAGIC,(b.stock('other_samples','Otheractualindexedlibraryallocations in sourcepool; sourcebodiesandamountsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.split('pool_loading','Actualsequencingloadingallocationfrompreparedmultisamplepool; amountunassigned','Selectedpooledlibraryloadingaliquot','Untransferredpreparedsequencingpool')
    lib=g.product('library','Published-source SRTRBsequencinglibrary; source spatial primers declared; full amplicons/index primers unassembled',K.LIBRARY);g.react('load','SourceNovaSeq6000SP; R135/I6/R2259',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(read_lengths={'R1':35,'I1':6,'R2':259},limitation='The published run lengths are reported here. Data S1 supplies the spatial primer/probe inventories; complete index-primer preparations, index orientation, occupied transcript sites and full amplicons remain unassigned. Run lengths alone do not establish molecular products.'))
    g.add_source_item(SourceItem('published_scope',AUTHORITY,'Final published article supports this selected SR-TRB preparation; published SM and separate Data S1 retrieved; Visium Rev A full-length cDNA and PN-2000089 primers are now inherited from the source-owned upstream declaration; historical capture and selected index-primer preparations remain unresolved. No whole-version equivalence or confirmed absence.'))
    g.add_source_item(SourceItem('analysis_metadata','Final published article, LR Spatial VDJ read processing, text1685; SR clonotype analysis1693–1700','Final article verifies16nt spatial barcode/12nt UMI in LR parsing. These are analysis coordinates, not whole primer sequences or a verified SR primer-binding geometry. SR processing groups by barcode plus UMI.'))
    g.add_source_item(SourceItem('preamp_programme',PREAMP_PROGRAMME,'Retainactualprintedvalues, executedpreamplificationcountunassigned. OtherselectedSRcyclecountsarelocallyexplicit.'))
    g.declare_reviewed_scope(MAIN_COMMON_TARGETS,
        'Final published article: target enrichment Methods1625–1642 and selected SR-TRB Methods1645–1668. Functional stock components are copresent classes, not additional doses or full formulation claims.')
    for replicate in SOURCE_REPLICATES:
        g.declare_reviewed_scope(tuple('rep'+str(replicate)+'_'+name for name in REPLICATE_MAIN_TARGETS),
            'Final SR-TRB Methods1645–1668: four distinct outer/inner replicates, analytical fractions, separate cleanup and10ng per replicate pooled at indexing. Unknown quantities remain unassigned.')
        g.declare_reviewed_scope(tuple('rep'+str(replicate)+'_'+name for name in DATA_S1_REPLICATE_TARGETS),
            'Published Data S1: PartRead1 worksheet and whole SR-Spatial VDJ outer/inner TRBV rows. Sequence-only cells do not establish terminal modifications or occupied transcript sites.')
    g.declare_reviewed_scope(DATA_S1_COMMON_TARGETS,
        'Published Data S1: Hyb. capture - Probes and Hyb. capture - Blocking Oligos. Actual source rows retained, including repeated sequences; pool molarity and end modifications unassigned.')
    for item, target in (('published_scope','source'),('analysis_metadata','source'),('preamp_programme','preamp')):
        g.cover_source_item(SourceCoverage(item,'modelled',target))
    g.annotate(SourceDiscrepancy(
        'preamp', 'Preamplification cycle count',
        (SourceClaim('Final article, target enrichment programme, item5 (text1632)',
                     'Repeat steps 2 to 5 six times'),
         SourceClaim('Same published item', 'for a total of five cycles')),
        'The preamplification cycle count remains unassigned.',
        'Both values are printed in the final published Methods; no corrected value is supplied by this model.',
        'The intended count is unresolved. The selected postcapture, outer, inner and indexing programmes are reported separately.',
        AnnotationSlot.BEFORE_WORKFLOW))
    g.annotate(SourceDiscrepancy(
        'rep1_outer_v','Outer TRBV worksheet label',
        (SourceClaim('Data S1 worksheet name and primer IDs','SR-Spatial VDJ outer TRBV; individual TRBV primer IDs'),
         SourceClaim('Same worksheet, interior heading A1','SR- Spatial VDJ Outer TRAV primers')),
        'Use the TRBV worksheet and individual row identities for the selected TRB preparation.',
        'The final Methods describe TRBV outer primers for this branch; the published worksheet bodies are preserved verbatim.',
        'The interior TRAV heading is inconsistent with the worksheet and row labels; it is not silently rewritten.',
        AnnotationSlot.BEFORE_WORKFLOW))
    return g.require_publishable()


def primer_example():
    # Separate source reagent orders, no invented transcript binding coordinates.
    sc=Scene()
    sc.strand('PartRead1',[Segment('PartRead1; also Visium cDNA forward',oligos.PARTREAD1,'r1')])
    sc.strand('Visium_cDNA_reverse',[Segment('PN-2000089 partial TSO reverse',visium_rev_a.REVERSE,'tso')])
    sc.strand('TRBV10-1_outer',[Segment('Outer V-specific primer',oligos.OUTER_SOURCE['TRBV10-1_outer'].upper())])
    name='TRBV10-1_inner';seq=oligos.INNER[name]
    sc.strand(name,[Segment('Nonbinding first-cycle 5′ Read2 handle',oligos.illumina.TRUSEQ_READ2,'r2'),
                    Segment('First-cycle 3′ V-specific binding part',seq[len(oligos.illumina.TRUSEQ_READ2):])])
    return sc


def partread1_binding():
    # A denatured strand from the actual upstream duplex; the other strand is
    # retained outside this binding-site view, not chemically deleted.
    sc=Scene()
    sc.strand('denatured_template',complement_segments(list(visium_rev_a.full_cdna())),rev=True,label='One strand of denatured Visium Rev A full-length cDNA')
    sc.anneal('PartRead1',[Segment('universal primer',oligos.PARTREAD1,'r1')],to='denatured_template',pair=('universal primer',"partial Read 1'"),label='PartRead1;3′ extension into tagged cDNA')
    sc.arrow('PartRead1','DNA synthesis into barcode/UMI/transcript side')
    return sc


def primer_design_rows():
    name='TRBV10-1_inner';whole=oligos.INNER[name];handle=oligos.illumina.TRUSEQ_READ2;target=oligos.INNER_SOURCE[name][1].upper()
    # Source columns are independently retained; this is a design explanation
    # of their existing checked assembly, not a transcript/amplicon assertion.
    if whole!=handle+target:raise ValueError('Inner primer does not match its source handle/target columns')
    code=lambda value:'<code>'+html.escape(value)+'</code>'
    return (
        ('cDNA amplification, before and after capture','Forward: '+code(visium_rev_a.FORWARD)+'<br>Reverse: '+code(visium_rev_a.REVERSE),
         'PN-2000089 forward binds the Read1 end; reverse binds the TSO-derived end. These retain full-length cDNA. Their whole sequences come from the cited Visium Rev A appendix.'),
        ('PartRead1, both receptor PCR rounds',code(oligos.PARTREAD1),
         'Use this whole '+str(len(oligos.PARTREAD1))+'-base order. It is the same sequence as the upstream cDNA forward primer, but pairs here with a V-specific primer. It contains no new sample index, spatial barcode or UMI.'),
        ('TRBV10-1 outer, receptor selection',code(oligos.OUTER_SOURCE['TRBV10-1_outer'].upper()),
         'The whole '+str(len(oligos.OUTER_SOURCE['TRBV10-1_outer']))+'-base source order is V-specific. No Read2 handle is prepended in this source outer primer. Pair with PartRead1.'),
        ('TRBV10-1 inner, handle assembly',code(handle)+' + '+code(target),
         str(len(handle))+'-base 5′ handle + '+str(len(target))+'-base 3′ V-specific portion = '+str(len(whole))+' bases. Order the concatenation in this direction. The 3′ portion binds the first template; the 5′ handle is copied into later products. Pair with the same PartRead1.'),
        ('TRBV10-1 inner, complete sequence',code(whole),
         'This is the published whole order, not its reverse complement. The shared handle is compatible with the source TruSeq indexing route; it does not specify which index-primer stock was used.'),
    )


def oligo_tables():
    groups=[]
    for title, rows in (
        ('Upstream cDNA primers and PartRead1',[("PN-2000089 forward / PartRead1",visium_rev_a.FORWARD),("PN-2000089 reverse",visium_rev_a.REVERSE)]),
        ('Outer TRBV primers',[(name,seq.upper())for name,seq in oligos.OUTER_SOURCE.items()]),
        ('Inner TRBV primers',list(oligos.INNER.items())),
        ('Capture probes',[(pool+' row '+str(i),seq)for i,(pool,seq)in enumerate(oligos.CAPTURE_SOURCE,1)]),
        ('Blocking oligos',list(oligos.BLOCKER_SOURCE.items())),
    ):
        groups.append('<details><summary>'+html.escape(title)+' ('+str(len(rows))+')</summary>'+
                      table(('Name','Sequence (5′→3′)','nt'),[(html.escape(name),'<code>'+html.escape(seq)+'</code>',str(len(seq)))for name,seq in rows])+'</details>')
    return ''.join(groups)


def receptor_class(stage):
    """Expected region order for a selected TRB product, never a literal locus.

    Primer-contributed ends are computed from actual source oligos. The variable
    transcript stays a role token: no allele, V-site spacing or full amplicon is
    manufactured. These inferred classes do not replace the graph's unresolved
    molecular transitions or establish source-specific PCR occupancy.
    """
    if stage not in {'captured','outer','inner','indexed'}:
        raise ValueError('unknown Spatial VDJ product class')
    prefix=list(visium_rev_a.full_cdna())[:6]
    # This variable region is explicitly not assigned a physical nucleotide
    # length. Its source-oriented class order follows 3′ capture of receptor RNA.
    body=[role('TRB: constant / J / CDR3 / retained V',inferred=True)]
    if stage=='captured':
        parts=prefix+body+list(visium_rev_a.full_cdna())[-2:]
    elif stage=='outer':
        parts=prefix+body+[Segment('outer-primer-derived V end',
            revcomp(oligos.OUTER_SOURCE['TRBV10-1_outer'].upper()),'r2',inferred=True)]
    else:
        target=oligos.INNER_SOURCE['TRBV10-1_inner'][1].upper()
        parts=prefix+body+[
            Segment('inner-primer-derived V end',revcomp(target),'r2',inferred=True),
            Segment('copied inner 5′ Read2 handle',revcomp(oligos.illumina.TRUSEQ_READ2),'r2',inferred=True)]
        if stage=='indexed':
            parts=[role('P5 / indexing-primer extension; exact order unassembled',inferred=True)]+parts+[
                role('sample index / P7 extension; exact order unassembled',inferred=True)]
    return Construct(parts,'INFERRED — '+stage+' TRB molecule class; not an occupied transcript reconstruction')


def step_scene(stage,*,primers=False):
    con=visium_rev_a.full_cdna() if stage=='input' else receptor_class(stage)
    bottom=complement_segments(list(con))
    sc=Scene.duplex(list(con),bottom=bottom,
        unpaired=tuple(seg.name for seg in bottom if seg.is_role_token()),
        label='Full-length cDNA' if stage=='input' else 'INFERRED — '+stage+' TRB class')
    sc.labels('top')
    if primers:
        # These are independent reagent orders. Do not place a V primer against
        # a made-up complementary template and claim actual binding validation.
        sc.strand('partread1',[Segment('PartRead1',oligos.PARTREAD1,'r1')],label='Incoming PartRead1 → barcode / UMI / receptor')
        if primers=='outer':
            sc.strand('v_primer',[Segment('TRBV10-1 outer',oligos.OUTER_SOURCE['TRBV10-1_outer'].upper(),'r2')],
                label='Incoming outer V primer → constant-region / tag side')
        elif primers=='inner':
            sc.strand('v_primer',[
                Segment('5′ Read2 overhang',oligos.illumina.TRUSEQ_READ2,'r2'),
                Segment('3′ V-specific binding part',oligos.INNER_SOURCE['TRBV10-1_inner'][1].upper())],
                label='Incoming inner V primer; 3′ target binds, 5′ handle overhangs')
        else:
            sc.strand('tso_primer',[Segment('PN-2000089 reverse',visium_rev_a.REVERSE,'tso')],
                label='Incoming partial TSO primer → full receptor / tag side')
        sc.note('v_primer' if primers in {'outer','inner'} else 'tso_primer',
            'Separate reagent rows; V-allele coordinates and occupied primer complexes remain unassembled.')
    return sc


def molecular_walkthrough():
    """One connected explanatory Workflow, complementary to the process audit.

    This shows source-supported transformations and conditional molecular
    classes. It is not a second validated preparation graph or a simulation.
    """
    w=Workflow(MolecularState('0 · Input: Visium CG000239 Rev A full-length cDNA, before fragmentation',
                              tuple(step_scene('input',primers='cdna').rows())),track_graph=False)
    def add(action,stage,name,explanation,primers=False):
        sc=step_scene(stage,primers=primers)
        sc.note('bottom',name)
        w.react(action,sc.rows(),name=name,note=explanation)
    add('1 · Preamplify separate 10 ng subsamples with the cDNA primer pair','input',
        'Retain both full-length cDNA handles and the original spatial barcode / UMI',
        'Take 10 µl of each Visium cDNA stock, split into 10 ng template subsamples, and amplify each separately. Each reaction: 25 µl 2× KAPA HiFi HotStart ReadyMix KK2602, 7.5 µl cDNA primers and 17.5 µl sample/water. Programme: 98°C 3 min; 98°C 15 s / 63°C 30 s / 72°C 2 min; 72°C 1 min final. The source contradicts itself on cycle count (“six times” / “five cycles”); the intended count is unresolved.')
    add('2 · Pool completed reactions, purify and concentrate','input',
        'Full-length cDNA retained; no fragmentation or new sequence added',
        'Pool only after the separate PCRs. Purify with 0.6× SPRIselect and elute each library in 40 µl Qiagen EB. The source then delegates concentration and capture setup to IDT xGen DNA-library capture version 4, optional Appendix A, then Plate protocol step 13. Those delegated operations are not replaced by a current kit recipe.')
    add('3 · Hybridize IG + TRA/TRB constant-region baits for 16 h; capture and wash','captured',
        'INFERRED — enriched full-length receptor class; the spatial tag remains connected',
        'Data S1 supplies the probe pools and universal-adapter blockers. Capture selects receptor molecules without intentionally trimming the ends. The selected experiment uses IG:TR 1:3: TR sequences were triplicated in the probe order and equal volumes of the produced IG/TR pools were mixed. This is not a threefold TR-pool volume. Exact occupied bait/template complexes and linked historical capture chemistry remain unassembled.',primers='cdna')
    add('4 · Postcapture PCR: 14 cycles with the same cDNA primer pair','captured',
        'INFERRED — amplified receptor-enriched full-length cDNA; clean before V-specific PCR',
        'This is the second experiment’s single 14-cycle branch, with the common cDNA-PCR reagents/programme, followed by AMPure purification. The first experiment’s 7-cycle then 5-cycle alternative is not combined into this route. The receptor and spatial tag still share a full-length molecule.',primers='outer')
    add('5 · In each of four replicates, outer TRBV + PartRead1 PCR, 12 cycles','outer',
        'INFERRED — outer TRB product: tag end retained, TSO-side end replaced by the V-primer boundary',
        'Use 2.5–5 ng captured cDNA per replicate. PartRead1 binds the universal tag-side handle and extends toward the receptor. The outer V primer binds upstream of CDR3 in transcript orientation and extends toward the constant region and tag. The selected product retains spatial barcode, UMI, CDR3 and constant-region sequence; it no longer needs the distant TSO end. The drawing shows a TRBV10-1 primer-derived end as one source example, not the full pool’s occupancy.')
    add('6 · Clean and quantify each outer product; transfer 5 ng to inner PCR','outer',
        'INFERRED — each replicate remains separate; no new DNA segment is introduced',
        'Source cleanup: 0.6× AMPure, two 80% ethanol washes, EB elution after 15 min at 37°C; Qubit and Bioanalyzer. Use qPCR to calibrate cycling. Transfer 5 ng from each replicate to its own inner reaction.',primers='inner')
    add('7 · Inner TRBV + the same PartRead1 PCR, 8 cycles','inner',
        'INFERRED — bounded inner-product class after copying both primer ends',
        'Only the V-side primer changes, so this is seminested PCR. Its 3′ V-specific portion binds closer to CDR3; its 5′ Read2-compatible tail initially overhangs the template. Extension retains that tail on the primer strand; subsequent template copying makes its complement. The duplex shown is the later bounded product class, not an immediate first-cycle duplex. The diagram retains the source handle in the correct strand polarity; actual V sites and allele-specific products remain unassembled.')
    add('8 · Clean and quantify; pool 10 ng from each of the four inner replicates','inner',
        'INFERRED — one sample-level pool; original spot barcodes and UMIs still distinguish molecules',
        'Use the source AMPure/ethanol cleanup and Qubit/Bioanalyzer checks. Pool 10 ng per replicate at indexing, not earlier. All four replicates receive one sample index and cannot be distinguished as technical replicates afterwards.')
    add('9 · TruSeq index PCR, 8 cycles; AMPure cleanup and QC','indexed',
        'INFERRED — final library architecture; unresolved whole index-primer arms are explicit',
        'Indexing adds the external sequencing/sample-index arms; it does not create new spatial barcodes or UMIs. The final Methods delegate the selected TruSeq index primers to reference 6. Their actual whole orders, index polarity and final primer landing remain unassembled; bracketed arms are not invented adapter sequences.')
    add('10 · Pool indexed samples and sequence on NovaSeq 6000 SP','indexed',
        'INFERRED — same library enters sequencing: Read1 35 nt / index 6 nt / Read2 259 nt',
        'Read1 recovers spatial-tag information; the long Read2 is designed to cover CDR3 and part of the constant region from the V side. The published run lengths are exact; this region-level class does not establish the first sequenced base or every SR cycle boundary. Map receptor calls to the retained Visium barcode and collapse by barcode plus UMI.')
    return w


def render_page():
    g=protocol_graph()
    paper='<a href="https://www.science.org/doi/full/10.1126/science.adf8486">Published article</a>'
    data='<a href="https://www.science.org/doi/suppl/10.1126/science.adf8486/suppl_file/science.adf8486_data_s1_to_s5.zip">Data S1 (within Data S1–S5)</a>'
    reads=table(('Read', 'Cycles', 'How to interpret it'), (
        ('Read 1', str(g.sequencing['library'].read_lengths['R1']),
         'Short spatial-tag read. The published run recipe is explicit; an exact cycle-by-cycle SR binding and base map is not assembled here.'),
        ('Index read', str(g.sequencing['library'].read_lengths['I1']),
         'Sample demultiplexing. The chosen index-primer orders and index orientation remain unassigned in this model.'),
        ('Read 2', str(g.sequencing['library'].read_lengths['R2']),
         'Long receptor read. The authors designed the asymmetric run to cover the entire CDR3 and part of the constant region from the 5′ end.'),
    ))
    identifiers=table(('Information', 'Meaning', 'Interpretation limit'), (
        ('Spatial barcode', 'Links a receptor molecule to a Visium tissue spot and the matched gene-expression map.',
         'A spot is a spatial sampling unit, not an individual-cell identifier. The final article’s LR parsing extracts 16 bases after the adapter anchor.'),
        ('UMI', 'Groups amplified reads derived from the same tagged molecule within a spot.',
         'The LR parsing extracts the next 12 bases. The SR pipeline groups by the combination of spatial barcode and UMI; UMI counts are not cell counts.'),
        ('CDR3 and receptor annotation', 'Identify and compare receptor clonotypes using the receptor sequence and MiXCR calls.',
         'A clonotype label is a computational assignment, not a sample index, antigen-specificity measurement or physical pairing of receptor chains.'),
        ('Sample index', 'Separates prepared tissue-section libraries multiplexed on the sequencer.',
         'All four TRB PCR replicates for an individual sample receive the same index.'),
    ))
    branches=table(('Route', 'Receptor coverage and purpose', 'Scope on this page'), (
        ('SR Spatial VDJ', 'TRA or TRB, using seminested V-region PCR and asymmetric Illumina reads.',
         'The TRB branch is represented: 12-cycle outer PCR, 8-cycle inner PCR. TRA uses a 14-cycle outer PCR and different V primers.'),
        ('LR Spatial VDJ', 'Full-length IG and TR transcripts, prepared as PacBio SMRTbell libraries; used for IG somatic hypermutation and class-switch analyses.',
         'Context only. Its adapter ligation and long-read preparation are not modelled by the SR workflow.'),
        ('Bulk / validation PCR', 'Separate experiments validate receptor sequences or evaluate capture.',
         'These are not interchangeable spatial libraries; the supplementary Fig. S1F validation PCR loses the spatial barcode.'),
    ))
    return assemble_page(g,[head('Spatial VDJ — spatial receptor sequencing, short-read TRB'),
        '<div class="wrap" data-protocol-complete="true"><h1>Spatial VDJ</h1>',
        '<p class="subtitle">Map lymphocyte receptor clones back to their tissue locations.</p>',
        info('<a href="https://doi.org/10.1126/science.adf8486">Engblom, Thrane, Lin et al., <i>Science</i> (2023)</a>. This schematic follows the published <strong>short-read TRB</strong> preparation, starting from full-length Visium cDNA before fragmentation.'),
        '<nav aria-label="On this page"><a href="#principle">What it measures</a> · <a href="#preparation">Preparation</a> · <a href="#primers">Primers</a> · <a href="#readout">Readout</a> · <a href="#limits">Scope and limits</a></nav>',
        '<h2 id="principle">Connect receptor identity to tissue position</h2>',
        '<p>Spatial VDJ adds immune-receptor sequencing to a spatial gene-expression experiment. A receptor’s complementarity-determining region 3 (<strong>CDR3</strong>) helps distinguish lymphocyte clones. The spatial barcode places that sequence on the tissue image, so clone distributions can be compared with gene-expression patterns and anatomical regions.</p>',
        '<p>The crucial starting point is <strong>full-length cDNA before fragmentation</strong>. The receptor sequence and the spatial tag are still connected at that stage. Constant-region hybridization capture enriches rare receptor transcripts; the short-read route then uses seminested PCR to bring the CDR3 region within reach of Illumina sequencing. '+paper+', Results “Spatial VDJ” and target-enrichment/SR Methods.</p>',
        HandoffSection(),
         '<h2 id="preparation">Step by step: from full-length cDNA to the TRB library</h2>',
        '<p>Follow the connected molecular states below. The top strand runs from the spatial-tag end toward the receptor’s V end; its receptor sequence therefore appears in the reverse order to the original RNA. A variable <strong>constant / J / CDR3 / V</strong> region replaces the unknown transcript, with no invented allele or measurable length. Each reagent row gives a source primer, not a fabricated occupied V-gene binding site.</p>',
        caveat('INFERRED — the walkthrough shows expected molecular classes, not validated full amplicons. The upstream Visium duplex and input primer orders are source declarations. The selected receptor products retain those source tag features; their V-site placement, amplicon sequence and index-primer arms remain unresolved. Dotted inferred sequence and bracketed regions preserve those limits. The separate process graph below records the actual partial-validation boundary.'),
        workflow_panel(molecular_walkthrough(),caption='Step-by-step SR-TRB preparation: source input, receptor capture, outer PCR, seminested handle addition and indexed library class'),
        table(('Receptor / index PCR','Published programme'),(
            ('Outer TRB, 12 cycles','98°C 5 min; 12 × (98°C 20 s, 65°C 30 s, 72°C 90 s); 72°C 7 min.'),
            ('Inner TRB, 8 cycles','98°C 5 min; 8 × (98°C 20 s, 72°C 30 s, 72°C 90 s); 72°C 7 min. Both consecutive 72°C entries are retained as printed.'),
            ('Indexing, 8 cycles','98°C 3 min; 8 × (98°C 20 s, 60°C 30 s, 72°C 90 s); 72°C 5 min.'))),
        '<p><strong>Why capture IG for a TRB experiment?</strong> The authors cocaptured IG and TR to provide enough product for the capture protocol despite low TR transcript abundance. They recommend both probe pools for either sequencing route. The selected 1:3 ratio comes from triplicating TR probe sequences in the order and mixing equal volumes of the resulting IG and TR pools; it does not mean adding three volumes of TR pool. '+paper+', “Target enrichment of spatial cDNA”.</p>',
        '<h2 id="primers">How to reconstruct the primers</h2>',
        '<p>All ordering sequences below are written <strong>5′→3′</strong>. The upstream Visium Rev A appendix supplies the cDNA primer pair; the receptor-specific orders come from published Data S1. Gene-specific target portions and shared 5′ tails are kept separate so that a handle is not mistaken for a genomic binding site.</p>',
        table(('Primer or stage','Sequence or construction,5′→3′','Design and binding role'),primer_design_rows()),
        panel(partread1_binding().rows(),caption='Computed PartRead1 binding on a denatured strand of the imported Visium Rev A cDNA. Barcode and UMI lie on the extension side and are retained in receptor-directed products. The other original strand is not shown in this site view; no TRBV coordinate or full receptor amplicon is asserted.'),
        '<p><strong>PartRead1</strong> recognizes the universal partial Read1 sequence of the spatial cDNA and is reused in both receptor-specific PCRs. <strong>Outer TRBV</strong> primers target V genes upstream of CDR3. <strong>Inner TRBV</strong> primers target closer to CDR3 and carry a shared Read2-compatible 5′ handle for subsequent TruSeq indexing. One matching V-gene example is shown below; each reaction uses the published primer pool, not just this pair. '+data+', SR-Spatial VDJ worksheets.</p>',
        '<p><strong>Can reads resolve the index-primer ends?</strong> Raw tonsil FASTQs are deposited under controlled access at <a href="https://ega-archive.org/datasets/EGAD00001011062">EGA</a>. The public <a href="https://zenodo.org/records/7961605">Zenodo archive</a> also contains processed TSV exports. We checked only their <code>targetSequences</code> fields in the first 100 TRB rows, finding no informative P5/P7/Read1/Read2 adapter read-through; this does not establish absence in the full dataset. Saved-input VDJCA contents and raw EGA reads were not inspected. Upstream trimming and UMI collapsing still limit interpretation. Even a copied terminal sequence would not uniquely establish the original primer synthesis order or modifications; the exact index-primer arms remain unresolved.</p>',
        panel(primer_example().rows(),caption='Independent source primer orders: the universal pair and one TRBV outer/inner example. Only the3′ target-specific part of an inner primer participates in its first-cycle receptor binding; its shared 5′ handle is carried as an overhang.'),
        '<p><strong>Reusing the design.</strong> Choose the matching outer/inner gene entry from the published pool: the same gene name identifies a pair, not identical binding positions. Outer targets lie 5′ of CDR3; inner targets lie closer to CDR3. Keep the source Read2 handle at the 5′ end of an inner primer, with its V-specific sequence at the3′ end. No universal sequence alone selects a receptor. Re-targeting requires the relevant transcript/allele and actual binding coordinates; the displayed variable transcript cannot establish those sites, amplicon length or specificity. Gene-specific annealing should be assessed from the binding portion rather than treating the nonbinding5′ handle as a first-cycle duplex.</p>',
        '<p><strong>What is inherited, and what is added?</strong> Spatial barcode and UMI come from the Visium starting cDNA. The outer/inner primers do not create new ones. The inner primer adds the Read2-compatible handle. The later index PCR adds sample demultiplexing information; its exact whole primer orders and index orientation still need reconstruction from the source-linked indexing method or verified read-derived evidence.</p>',
        '<p><strong>Capture probes</strong> target IG and TRA/TRB constant regions. <strong>Blocking oligos</strong> occupy universal adapter sequences to reduce off-target capture. The Read1+polyT blocker contains a degenerate N window and a polyT tract; those N bases are a reagent pattern, not a declared barcode or UMI for a particular molecule.</p>',
        panel(oligos.pool_scene('blockers').rows(),caption='The two published capture-blocker patterns from Data S1. N and V retain their source degeneracy; terminal modifications are not assigned from sequence-only cells.'),
        '<details><summary>All published oligo sequences</summary><p>Whole source orders are retained, including repeated source rows. Sequences alone do not establish terminal modifications, occupied transcript binding sites or a complete library construct.</p>',oligo_tables(),'</details>',
        '<h2 id="readout">What the reads and identifiers tell you</h2>',reads,
        '<p>The final article gives an explicit 16-base spatial barcode followed by a 12-base UMI in its <em>long-read processing</em> description. That verifies the tag interpretation, but does not by itself establish every short-read cycle boundary or a complete Visium primer structure. '+paper+', “LR Spatial VDJ read processing” and “SR Spatial VDJ clonotype analysis”.</p>',identifiers,
        '<p><strong>From reads to a tissue map.</strong> The published SR analysis filters low-quality reads, requires a valid Visium barcode and the inner-V primer sequence, and collapses reads by barcode plus UMI. For tonsil repertoire analysis, each UMI group must have at least three supporting reads. MiXCR supplies CDR3/clonotype calls, which become a clone-by-spatial-barcode count matrix and are matched to the gene-expression map. The authors also filter index-hopping candidates and a small number of suspected overamplification artifacts. These analysis choices affect the reported clone map; sequencing yield alone is not a count of clones or cells. '+paper+', “SR Spatial VDJ clonotype analysis”.</p>',
        '<h2 id="limits">Branches and interpretation limits</h2>',branches,
        '<p><strong>Spatial colocalization does not physically pair receptor chains.</strong> Spatial VDJ captures individual receptor transcripts. The paper predicts IG heavy/light pairing computationally with “repair” and evaluates predictions using matched single-cell VDJ data. A TRB map alone does not identify its paired TRA chain or establish antigen specificity. '+paper+', receptor-pairing Results.</p>',
        caveat('This page represents the published preparation order and source oligos. The pre-fragmentation starting duplex and cDNA primer pair now come from the exact cited Visium Rev A declaration. Source-specific index-primer orders, terminal modifications of the Data S1 reagents, occupied receptor sites, full amplicons and exact final sequencing-primer landing remain unassembled. These are model boundaries, not a claim that the publication omits those details.'),
        '<details><summary>Read-plan detail</summary>',PlanLimitationSection('library'),'</details>',
        '<h2>Published source discrepancies</h2>',AnnotationSlot.BEFORE_WORKFLOW,
        '<details><summary>Full reaction workflow, physical fractions and validation</summary><p>The expanded workflow retains separate PCR replicates, cleanup fractions, analytical aliquots and stock remainders. Its unresolved transitions expose where molecular products or reagent geometry have not been assembled.</p>',
        WorkflowSection('Source-qualified SR-TRB reaction workflow'),'</details>',
        '<p class="research-notes"><a href="01_spatial_vdj.html">Research notes and source-reading scope</a></p>',
        '</div>'])


def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    # The explanatory product class must preserve upstream identifier identity
    # and the polarity of a handle copied from the incoming inner primer.
    upstream_tags={seg.name:seg.feature for seg in visium_rev_a.full_cdna() if seg.feature}
    for stage in ('captured','outer','inner','indexed'):
        model=receptor_class(stage)
        c(stage+' class retains the same upstream tag features',
          {seg.name:seg.feature for seg in model if seg.feature},upstream_tags)
        c(stage+' unknown receptor region is not assigned an invented nucleotide length',
          next(seg for seg in model if seg.is_role_token() and seg.name.startswith('TRB:')).length_bp,None)
    model=receptor_class('inner')
    c('Inner product carries the full source primer on the opposite strand',
      revcomp(model.top()[-len(oligos.INNER['TRBV10-1_inner']):]),oligos.INNER['TRBV10-1_inner'])
    c('Outer-specific primer sequence is not silently inherited at the inner boundary',
      list(model)[-2].top!=revcomp(oligos.OUTER_SOURCE['TRBV10-1_outer'].upper()),True)
    upstream=g.imports['source']
    c('Starting cDNA is imported from the exact source Visium revision',upstream.reference.version,'CG000239 Rev A')
    c('Starting cDNA is the upstream produced cleanup output',upstream.reference.producer,'cdna_cleanup')
    c('Starting cDNA cannot be a fragmented sequencing library',g.materials['source'].chemical_state['fragmented'],False)
    c('Both source-named cDNA primer stocks are now assembled',all(len(g.materials[key].state.strands)==2 for key in ('preamp_cdna_primers','postcapture_cdna_primers')),True)
    c('Four actual outerTRB template allocations plus stockremainder',len(g.steps['outer_allocate'].outputs),5)
    c('Physicalfourreplicatepoolonlyafterselectedinnerproducts',tuple(p.material for p in g.steps['pool_replicates'].inputs),tuple('rep'+str(i)+'_select10ng_selected'for i in range(1,5)))
    c('OneindexPCRconsumesactualreplicatepool',g.steps['index'].inputs[0].material,'index_pool')
    c('Postcaptureuses14cyclebranchaftercohybridization',g.steps['postcapture'].inputs[0].material,'capture_product')
    c('Eachinnerusesitsown5ngouterallocation',all(g.steps['rep'+str(i)+'_inner'].inputs[0].material=='rep'+str(i)+'_select5ng_selected'for i in range(1,5)),True)
    c('SelectedoutertemplatecomesafterqPCRanalyticalfractionremoved',g.steps['outer_allocate'].inputs[0].material,'outer_calibration_stock')
    c('InnerQCcannotentertheindexpool',all(g.terminals['rep'+str(i)+'_inner_qc_analytical'].disposition==Disposition.QC for i in range(1,5)),True)
    c('SamephysicalpreampReadyMixsuppliesknownenzymeandNT',tuple(g.steps['preamp_components'].outputs),('preamp_enzyme','preamp_nt','preamp_buffer'))
    c('ReadyMixcomponentproducerisoneaggregateinventoryinput',g.steps['preamp_components'].inputs[0].material,'preamp_kapastock')
    c('Variablepreampaggregateamountnotasingle25ulreaction','aggregate totalamountunassigned' in g.materials['preamp_kapastock'].chemical_state['volume'],True)
    c('IndexPCRdoesnotinheritKAPAidentity', 'not assigned' in g.materials['index_kapastock'].chemical_state['physical_stock'],True)
    c('Onlyselectedpooledloadingaliquotenterssequencer',g.steps['load'].inputs[0].material,'pool_loading_selected')
    c('Preparedsequencingpoolremainderiskeptuntransferred',g.terminals['pool_loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('Published IG pool preserves every source row',len(g.materials['ig_probes'].state.strands),174)
    c('Published TR pool preserves every source row',len(g.materials['tr_probes'].state.strands),20)
    c('All four actual outer allocations contain the source pool',all(len(g.materials['rep'+str(i)+'_outer_v'].state.strands)==len(oligos.OUTER_SOURCE)for i in range(1,5)),True)
    c('All four actual inner allocations contain the source pool',all(len(g.materials['rep'+str(i)+'_inner_v'].state.strands)==len(oligos.INNER_SOURCE)for i in range(1,5)),True)
    c('Both blockers retained with degenerate pattern',len(g.materials['capture_blockers'].state.strands),2)
    c('Known primer pools retain unavailable transcript geometry',all(g.steps['rep'+str(i)+'_'+stage].geometry.availability=='unavailable'for i in range(1,5)for stage in ('outer','inner')),True)
    c('Wholesourceconstructisnotfabricated',g.materials['library'].state,None)
    c('Readlengthsaresource35and259',g.sequencing['library'].read_lengths,{'R1':35,'I1':6,'R2':259})
    p=sources_root()/Path(__file__).resolve().parents[1].name/'science_adf8486_final_browser.html.txt'
    if have(c,Source(p,'https://www.science.org/doi/full/10.1126/science.adf8486')):
        text=p.read_text()
        fixtures=[
            ('Published preamplification ambiguity retained','Repeat steps 2 to 5 six times for a total of five cycles'),
            ('Equal pool volumes establish probe ratio','adding equivalent volumes of the produced pools'),
            ('Selected fourteen-cycle postcapture branch','reaction as described above, but with 14 cycles.'),
            ('Four actual replicates','Four replicate reactions were prepared per sample'),
            ('Inner PCR retains printed temperature','72°C: 30 s'),
            ('Ten ng per replicate enters indexing','10 ng of input material was used from each replicate'),
            ('Published run lengths','read 1, 35 nt; index read, 6 nt; read 2, 259 nt'),
            ('Primer sequences delegated to Data S1','see data S1 for primer sequences'),
            ('Final published spatial tag parsing','first 16 bases following the Truseq adapter for the spatial barcode and the 12 bases following that as the UMI'),
            ('Final published SR collapsing retains spot and molecule identity','both the UMI and spatial barcode sequence as the “UMI” region'),
            ('Final published receptor pairing limit','Spatial VDJ captures individual antigen receptor transcripts and cannot therefore define receptor pairs by itself'),
        ]
        for label,needle in fixtures:c(label,needle in text,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'science.adf8486_data_s1.xlsx.txt'
    if have(c,Source(p,'https://www.science.org/doi/suppl/10.1126/science.adf8486/suppl_file/science.adf8486_data_s1_to_s5.zip')):
        source_text=p.read_text().upper()
        c('Published complete PartRead1 transcription',oligos.PARTREAD1_SOURCE in source_text,True)
        c('All whole outer pool orders occur verbatim',all((name+'\t'+seq).upper() in source_text for name,seq in oligos.OUTER_SOURCE.items()),True)
        c('Independent whole inner orders match declared canonical assemblies',all(oligos.INNER[k]==whole.upper() for k,(whole,_)in oligos.INNER_SOURCE.items()),True)
        c('All whole inner source orders occur verbatim',all((name+'\t'+whole+'\t'+target).upper() in source_text for name,(whole,target)in oligos.INNER_SOURCE.items()),True)
        c('All capture probe rows occur verbatim',all((pool+'\t'+seq).upper() in source_text for pool,seq in oligos.CAPTURE_SOURCE),True)
        c('Both complete source blockers occur verbatim',all((name+'\t'+seq).upper() in source_text for name,seq in oligos.BLOCKER_SOURCE.items()),True)
    c.report()
if __name__=='__main__':selftest()
