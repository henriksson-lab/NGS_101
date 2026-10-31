from pathlib import Path
import html
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from chemdraw import Scene,Segment,feature
from pcr import PCRGeometry
import oligos_spatial_vdj as oligos
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition,AnnotationSlot,SourceClaim,SourceDiscrepancy
from page import head, info, caveat, table
from protocol_views import assemble_page, WorkflowSection, PlanLimitationSection
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
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True,chemical_state={})
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

AUTHORITY='Final published Science full article, DOI10.1126/science.adf8486, authorized browser HTML science_adf8486_final_browser.html.txt: Results1580/Methods1623–1669; supplementary index1767–1782. Earlier accepted-analysis-only16/12 metadata labelled locally. Published SM PDF title/Methods1–72 and source DataS1 spatial worksheets inspected after user supply; full kit/template geometry unresolved. Earlier accepted-only review is historical, not whole-version equivalence.'
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
        b.op(key,label,S.MAGIC,supplies,reason='Published amplification order; complete cDNA/index primer preparation and occupied products unassembled. Same actual ReadyMix components exposed; index supplies not assumed to be that premix.')

def cleanup(b,key,label):b.clean(key,label,extra=(b.stock(key+'_beads','Actualsource '+label+' bead allocation',K.SELECTION,R.SELECTION),))
def qc(b,key,label):
    g=b.g;y=g.product(key+'_analytical',label,K.QC);r=g.product(key+'_stock','Actual preparation remaining after analytical allocation; amountunassigned',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(y.id,Disposition.QC,'Actualanalyticalfraction excludedfromproduction'));b.x=r

def protocol_graph():
    b=B('Spatial VDJ published SR-TRB four-replicate partial preparation',AUTHORITY,'Actual full-length spatially barcoded Visium cDNA before fragmentation from published-source tonsil preparation',K.MOLECULAR);g=b.g
    g.materials['source'].chemical_state.update({'source_boundary':'Published Methods useCG000239RevA Visiumpreparation; its oligos/RT/TSOrecipe NOTimported.','spatial_barcode':'Earlier accepted-manuscript LR analysis843–865 identifies16nt after anchor; not independently final-verified or a complete oligo','UMI':'Earlier accepted-manuscript LR analysis843–865 identifies12nt after spatial barcode; not independently final-verified; origin/wholeconstructunassembled','scope':'StartsactualpreparedcDNA, notcompletefirst-strandchemistryorwholefinalconstruct'})
    b.split('select10ul','Actual10ulfromVisiumfull-lengthcDNAstock','Actualselected10ulcDNA','UntransferredVisiumcDNAstock')
    collection=g.product('preamp_collection','Actual collection of separate10ngtemplate PCR subsamples; number depends on measured concentration, not a physical pool',K.MOLECULAR);g.react('preamp_allocate','Splitselected10ulinto10ngtemplates, multiplicityunassigned',S.GENERIC,((b.x,R.TEMPLATE),),(collection,),changes_molecule=False);b.x=collection
    amp(b,'preamp','Commonpreamplification of separate10ngsubsamples, notpooledbeforePCR',(primer(b,'preamp_cdna_primers','Aggregateactual10xcDNAprimerallocations7.5ulPERPCR acrossvariablearray; number/totalamount/constituentwholeordersunassembled'),),PREAMP_PROGRAMME+' Source25ulReadyMix+7.5ulprimers+17.5ulsample/water;98C3min,98C15s/63C30s/72C2min;72C1minfinal.')
    b.op('preamp_pool','Source physical pooling of completedpreamplification subsamples',S.GENERIC,product='ActualpooledpreamplifiedspatialcDNA; parallelreactionmultiplicity/wholeproductsunassembled')
    cleanup(b,'preamp_clean','0.6xSPRIselect; elute40ulQiagenEB')
    b.op('concentrate','SourceIDTxGenv4optionalAppendixAAMPureconcentrationthenPlateProtocolstep13',S.MAGIC,(b.stock('capture_concentration','ActualAMPureconcentrationpreparation; delegatedamounts/completeconditionsunassembled',K.SELECTION,R.SELECTION),))
    b.op('capture','Source16hovernightIGandTRcohybridization/capture; selectedIG:TR1:3',S.MAGIC,(primer(b,'ig_probes','ActualIGconstantregionDiscoverypool; DataS1 whole pool now declared; end modifications unassembled'),primer(b,'tr_probes','ActualTRA/TRBconstantregionDiscoverypool; TRG/TRDnotincluded, DataS1 whole pool now declared; end modifications unassembled'),primer(b,'capture_blockers','Actualcustomuniversal-adapterblockingoligos; DataS1 whole pool now declared; end modifications unassembled'),b.stock('xgen_capture','ActualIDTxGenHybridizationandWashKit1080584 v4capture/washpreparation; selection/packagingcomponentsunassembled',K.SELECTION,R.SELECTION),b.cond('capture_ratio','TRbaitsequencestriplicatedinorder; equalvolumesofIG/TRproducedpoolsachieve1:3, NOTthreefoldTRpoolvolume.')),
         reason='Published-source hybridization capture established; baitmodifications/selectionligands/wholeboundRNA-derivedDNAgeometry and linkedIDTv4recipeunassembled. No assumed biotin/streptavidin or fresh vendor-default components.')
    amp(b,'postcapture','Selectedsecondexperiment:one14cyclepostcapturePCR',(primer(b,'postcapture_cdna_primers','Actual7.5ul10xcDNAprimerpreparation perPCR; wholeordersunassembled'),),'Source14cyclesbranch; firstexperiment7+5cyclebranchnotmerged; source25ulKAPAReadyMix/17.5ulsampleandwater, programmeotherwiseasreported.')
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
    g.sequence(lib,SequencingPlan(read_lengths={'R1':35,'I1':6,'R2':259},limitation='The published run lengths are reported here. Data S1 supplies the spatial primer/probe inventories; complete cDNA/index primer preparations, index orientation, occupied transcript sites and full amplicons remain unassigned. Run lengths alone do not establish molecular products.'))
    g.add_source_item(SourceItem('published_scope',AUTHORITY,'Final published article supports this selected SR-TRB preparation; published SM and separate Data S1 retrieved; historical vendor detail/whole cDNA and index primer preparations remain unresolved. No whole-version equivalence or confirmed absence.'))
    g.add_source_item(SourceItem('historical_analysis_metadata','Accepted author version Zenodo13148082, LR analysis843–865 only','16nt spatial barcode/12nt UMI are analysis metadata from earlier inspected accepted original, not final-verified whole primer sequences or binding evidence.'))
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
    for item, target in (('published_scope','source'),('historical_analysis_metadata','source'),('preamp_programme','preamp')):
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
    # Separate reagent orders: this view asserts no annealing or transcript site.
    sc=Scene()
    sc.strand('PartRead1',[Segment('PartRead1',oligos.PARTREAD1,'r1')])
    name='TRBV10-1_inner';seq=oligos.INNER[name]
    sc.strand(name,[Segment('Published Read2-compatible handle',oligos.illumina.TRUSEQ_READ2,'r2'),
                    Segment('Source TRBV10-1 target-specific portion',seq[len(oligos.illumina.TRUSEQ_READ2):])])
    return sc


def oligo_tables():
    groups=[]
    for title, rows in (
        ('Outer TRBV primers',[(name,seq.upper())for name,seq in oligos.OUTER_SOURCE.items()]),
        ('Inner TRBV primers',list(oligos.INNER.items())),
        ('Capture probes',[(pool+' row '+str(i),seq)for i,(pool,seq)in enumerate(oligos.CAPTURE_SOURCE,1)]),
        ('Blocking oligos',list(oligos.BLOCKER_SOURCE.items())),
    ):
        groups.append('<details><summary>'+html.escape(title)+' ('+str(len(rows))+')</summary>'+
                      table(('Name','Sequence (5′→3′)','nt'),[(html.escape(name),'<code>'+html.escape(seq)+'</code>',str(len(seq)))for name,seq in rows])+'</details>')
    return ''.join(groups)


def render_page():
    g=protocol_graph()
    overview=table(('Stage','Source-supported preparation'),(
        ('Starting material','Full-length spatially barcoded Visium cDNA, before fragmentation.'),
        ('Preamplification and capture','Amplify 10 ng subsamples separately, pool completed products, then capture IG/TR constant-region transcripts together overnight. The selected IG:TR probe ratio is1:3.'),
        ('Nested TRB amplification',str(len(SOURCE_REPLICATES))+' separate replicates: outer TRBV + PartRead1, followed by inner TRBV + PartRead1, with cleanup and QC between stages.'),
        ('Indexing','Pool 10 ng from each inner product at indexing; all four replicates receive the same sample index.'),
        ('Sequencing','NovaSeq 6000 SP; short R1 and long R2 recover the spatial tag and receptor sequence. Exact occupied primer geometry remains unassembled.'),
    ))
    reads=table(('Read','Cycles','Source scope'),[(html.escape(role),str(cycles),
        'Published run length; whole index order/orientation unassigned' if role=='I1' else
        'Published run length; exact read-base mapping not asserted')for role,cycles in g.sequencing['library'].read_lengths.items()])
    return assemble_page(g,[head('Spatial VDJ — short-read TRB'),
        '<div class="wrap" data-protocol-complete="true"><h1>Spatial VDJ</h1>',
        '<p class="research-notes"><a href="01_spatial_vdj.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1126/science.adf8486">Engblom et al., <i>Science</i> (2023)</a>. Spatially resolved immune-receptor sequencing. This page follows the published short-read TRB branch; the long-read IG/TR and separate validation preparations are outside its scope.'),
        caveat('The workflow and published oligo inventories are represented, with unresolved chemistry retained. Complete cDNA/index-primer preparations, terminal modifications and occupied transcript/product geometry remain unassembled. A linked page does not mean every reaction is molecularly validated.'),
        '<h2>Preparation</h2>',overview,
        '<h2>Published primer sequences</h2>',
        info('PartRead1 and one inner TRBV order are shown below as separate reagents. The full outer/inner pools, capture probes and blocking patterns are available in the tables; repeated source rows are retained.'),
        panel(primer_example().rows(),caption='Data S1 oligo orders; separate rows, not an annealed complex or a complete final library.'),
        oligo_tables(),
        '<h2>Read setup</h2>',reads,PlanLimitationSection('library'),
        AnnotationSlot.BEFORE_WORKFLOW,
        '<details><summary>Detailed reaction model and validation</summary>',
        WorkflowSection('Source-qualified SR-TRB reaction workflow'),'</details>',
        '</div>'])


def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
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
