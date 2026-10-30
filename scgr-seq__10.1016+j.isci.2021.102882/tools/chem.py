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


from dataclasses import dataclass,replace
import illumina as il
from protocol_graph import Terminal,Disposition
TITLE='scGR-seq — released lectin barcodes and same-cell RNA'
SOURCE='Published PMC8349903 Results88/Fig1 and Methods330–350;publishedTableS1rPVLrow andTableS15indexprimerfamilies'
LIMIT='Private partial source model. Photocleavable-linker chemistry is named in Results while Methods abbreviates DBCO-NHS;exact linker structure/released5prime remnant remain unassembled. Physical tag supernatant andcellpellet are distinct;noRNA/DNA idealpurity. HistoricalGenNextRamDAkitandcompletePCR/readgeometryunassembled.'
WHOLE_TAG='CGACGCTCTTCCGATCTCTGTCGCCCTGAACAACGGCAATGGCAGTCAGATCGGAAGAGCACAC'
LECTIN_ID='TCGCCCTGAACAACGGCAATGG'
@dataclass(frozen=True)
class LectinReporter:
    stage:str='azide'
    def __post_init__(self):
        if self.stage not in ('azide','conjugated','released'):raise ValueError('source reporter stage required')
    def scene(self):
        start=WHOLE_TAG.index(LECTIN_ID)
        sc=Scene();sc.strand('reporter',[Segment('source5primePCRarm',WHOLE_TAG[:start]),Segment('rPVLlectinidentity',LECTIN_ID,feature=feature('scgr.rpvl','feature_barcode','fixed')),Segment('source3primePCRarm',WHOLE_TAG[start+len(LECTIN_ID):])],label='Published rPVL DNA-tag example',mod5={'azide':'5prime azide','conjugated':'lectin/linker/triazole attachment;geometryunassembled','released':'photorelease remnant unassembled'}[self.stage])
        sc.note('reporter','Source '+self.stage+' reporter;not a glycan polymer,UMI or newly assignedcellbarcode. Fullorganiclinker/endformula unassembled.');sc.lectin_reporter=self;return sc

def index_primer(which):
    if which=='i5':return Construct([Segment('P5',il.P5),Segment('sourcecellwellindex','N'*8,placeholder=True,feature=feature('scgr.well_i5','cell_barcode','whitelist',whitelist='PublishedTableS15i5members',note='Sequencerindexalsoidentifiessinglecellwell;notUMI.')),Segment('sourcePCRarm',il.TRUSEQ_READ1+'CTG')])
    if which=='i7':return Construct([Segment('P7',il.P7),Segment('sourcecellwellindex','N'*8,placeholder=True,feature=feature('scgr.well_i7','cell_barcode','whitelist',whitelist='PublishedTableS15i7members',note='Sequencerindexalsoidentifiessinglecellwell;notUMI.')),Segment('sourcePCRarm','GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGA')])
    raise ValueError(which)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    lectin=g.input('lectin','Source lectin panel;representativerPVL,otherlectins/negativeIgGcontrolsnotcollapsedasonebinder',K.MOLECULAR)
    linker=g.input('linker','Source namedphotocleavableDBCO-NHSResults;MethodsabbreviatesDBCO-NHS,exactstructureunassembled',K.MOLECULAR)
    x=g.product('activated','DBCO-functionalizedlectin plusunreactedlinker/reagents',K.MOLECULAR)
    g.react('activate','Lectin+10xmolarNHSlinker20C1h dark',S.GENERIC,((lectin,R.TEMPLATE),(linker,R.STRUCTURAL)),(x,))
    x=op('quench',x,'Add1MTris10uL20C15mindark toquenchexcessNHS','SourceTrisquenchingcondition')
    x=clean('desalt',x,'SephadexG25 removeexcesslinker;noquantitativepurityassertion')
    dna=g.input('azide_tag','Source5prime-azideDNA;whole rPVLtagfromTableS1',K.OLIGO,state=LectinReporter().scene())
    out=g.product('conjugated','Lectin–linker–DNA conjugate classesplusfreeDNA/reactants',K.MOLECULAR,state=LectinReporter('conjugated').scene())
    g.react('click','Source4Ccouplingwith10xmolar5prime-azideDNA;reactiontimeunassembled',S.GENERIC,((x,R.TEMPLATE),(dna,R.STRUCTURAL)),(out,))
    x=clean('affinity',out,'Appropriatesugar–Sepharoseaffinitychromatography;removeunboundDNA/recoveractivelectinconjugates')
    x=clean('dialysis',x,'SourceTube-O-DIALYZERdialysis then10kDaMWCOconcentration;actualseparations')
    cells=g.input('cells','Source cellsinsuspension;singlecellpreparationmodel',K.BIOLOGICAL,primary=True)
    stained=g.product('stained','Cellswithboundlectinconjugates plusunboundpanel/buffer',K.MOLECULAR,state=LectinReporter('conjugated').scene())
    g.react('stain','1x10^5cells/41proteins0.5ug/ml4C1hPBS1%BSA',S.GENERIC,((cells,R.TEMPLATE),(x,R.STRUCTURAL)),(stained,))
    x=clean('cell_washes',stained,'Three1mLPBS/BSAwashes;retaincells')
    x=op('pick',x,'Countcells;TOPickvisualsingle-cellpickintoPCRtube;notbulk1x10^4branch')
    x=op('uv',x,'365nm15W15minphotorelease;cellandDNA/remnantsremaincopresent','SourceUVcondition',state=LectinReporter('released').scene())
    sup=g.product('barcode_supernatant','RecoveredreleasedDNA-barcodesupernatant;solublecarryoverpossible',K.MOLECULAR,state=LectinReporter('released').scene())
    pellet=g.product('cell_pellet','SamecellpelletretainedforRNA;lectin/linker/cellremnantsmayremain',K.MOLECULAR)
    spin=g.input('spin_condition','Sourcecentrifugation/supernatanttransfer;fullhandlingdetailsunreviewed',K.CONDITION)
    g.react('separate','Physicaltag-supernatant/cell-pelletseparationbeforeindependentlibraries',S.CLEANUP,((x,R.TEMPLATE),(spin,R.CONDITION)),(sup,pellet))
    f=g.input('i5','Publishedi5indexprimerfamily',K.OLIGO,state=index_primer('i5'))
    r=g.input('i7','Publishedi7indexprimerfamily',K.OLIGO,state=index_primer('i7'))
    enzyme=g.input('q5','NEBNextUltraIIQ5',K.CATALYST)
    pcr=g.product('barcode_pcr','IndexedDNA-tagampliconsandPCRmixture;organic5primeremnantnotassignedtotagcopies',K.MOLECULAR)
    g.react('barcode_pcr_step','Source20cycle98/65PCRwithcellwelli5/i7indexes;noRNAreverse-transcription',S.GENERIC,((sup,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(pcr,))
    lib=clean('barcode_cleanup',pcr,'SourceAMPureXPpurification,MultiNAassessment;ratiohistoricalvendorprotocolunreviewed')
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':26,'Read 2':26},limitation='SourceMiSeqv2 50cyclekit and26bpPEbothprinted;no repairto25bp,exactindexacquisition/primergeometryunassembled.'))
    rna=op('rna_prepare',pellet,'Sourcecell-pelletRNApreparation andGenNextRamDA-seqSingleCellKitcontinuation','HistoricalGenNextRamDAkit;exactextraction/primers/operationsunreviewed')
    rna=op('rna_qc',rna,'SourcePowerUpSYBR/CFXreal-timePCRlibraryquantificationwithDNAstandards')
    g.sequence(rna,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':151,'Read 2':151},limitation='SourceNovaSeq6000S4,151bpPE;wholehistoricalRamDA/librarypolymersandsequencingprimerlandingunassembled. NoDNAbarcodetransferintotheRNAbranchasserted.'))
    return incomplete(g)
def render_page():
    g=protocol_graph();return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_scgr.html">Research notes</a></p>'+''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ('barcode_cleanup_kept','rna_qc_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('ConjugatedDNAnolongerfree5azide',g.materials['conjugated'].state.lectin_reporter.stage,'conjugated')
    c('Releasedtagdoesnotclaimfresh5phosphate',g.materials['barcode_supernatant'].state.strands['reporter'].mod5,'photorelease remnant unassembled')
    c('TagPCRusesphysicallyrecoveredsupernatant',g.steps['barcode_pcr_step'].inputs[0].material,'barcode_supernatant')
    c('RNAcontinuesfromcellpellet',g.steps['rna_prepare'].inputs[0].material,'cell_pellet')
    c('NoRTinsertedinreleasedDNAtaglibrary',g.steps['barcode_pcr_step'].kind,S.GENERIC)
    c('Modelprivate',g.publishable,False)
    folder=sources_root()/Path(__file__).resolve().parents[1].name
    table=folder/'j.isci.2021.102882_supp_mmc2_TableS1.xlsx.txt'
    if table.exists():c('OriginalwholeDNAandlectinidentity',WHOLE_TAG+'\t'+LECTIN_ID in table.read_text(),True)
    else:c.skip('TableS1transcription','Externaloriginalabsent')
    table=folder/'j.isci.2021.102882_supp_mmc2_TableS15.xlsx.txt'
    if table.exists():
        text=table.read_text().replace('X','N');c('Bothoriginalindexprimerfamilies',all(x.top() in text for x in (index_primer('i5'),index_primer('i7'))),True)
    else:c.skip('TableS15transcription','Externaloriginalabsent')
    c.report()
