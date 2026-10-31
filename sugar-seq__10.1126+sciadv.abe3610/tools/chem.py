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
    g.declare_reviewed_scope(STAINING_SCOPE,SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

from dataclasses import dataclass
TITLE='SUGAR-seq — lectin-linked glycan reporters with RNA, TCR and epitopes'
SOURCE='Published PMC7895430 Results78–90/Fig1A and Methods160–162;main linked SI identityDC1 is distinct from downloaded figure-only zip'
LIMIT='Private source-limited operational graph. Glycans are detected via affinity-bound biotin-L-Pha/anti-biotin DNA reporter,not converted into DNA or covalently joined to RNA. Exact reporter/barcoding/library polymers,historical5prime kitversion and linked chemistrySI remain unassembled.'
@dataclass(frozen=True)
class ReporterComplex:
    lectin_bound:bool=False
    reporter_bound:bool=False
    def rows(self):
        sc=Scene();sc.strand('reporter',[Segment('Source anti-biotin antibody reporter','[DNA reporter polymer unassembled]',placeholder=True,feature=feature('sugar.glycan_reporter','feature_barcode','unknown'))],label='DNA tag belongs to anti-biotin antibody;not glycan backbone')
        sc.note('reporter','Source affinity complex:cellglycan—L-Pha carryingbiotin—anti-biotin antibody withDNA tag. Binding geometry/stoichiometry/recovery unassembled;not covalentglycan→DNA conversion. State lectinbound='+str(self.lectin_bound)+', reporterbound='+str(self.reporter_bound));return sc.rows()

STAINING_SCOPE=('cells','biotin_lpha','native_lpha','lectin_bind','lectin_wash_selection','lectin_wash','antibodies','conjugate_inventory','antibody_stain','antibody_wash_selection','antibody_wash')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    cells=g.input('cells','Source FACS-sorted tumor-infiltrating Tcells',K.BIOLOGICAL,primary=True)
    biotin=g.input('biotin_lpha','Source biotinylatedL-Pha carbohydrate-binding protein;biotinlocation/stoichiometry unassembled',K.MOLECULAR)
    native=g.input('native_lpha','Source nonbiotinylatedL-Pha,1:5mixture withlabelledL-Pha',K.MOLECULAR)
    x=g.product('lectin_stained','Cells withglycan-boundlectin plusunboundlectin;affinitybinding notglycanDNAconversion',K.MOLECULAR,state=ReporterComplex(True,False))
    g.react('lectin_bind','Source1ug/mlL-Pha mixture;20minonice',S.BINDING,((cells,R.TEMPLATE),(biotin,R.STRUCTURAL),(native,R.STRUCTURAL)),(x,),changes_molecule=False)
    x=clean('lectin_wash',x,'ThreecoldPBSwashes;retaincells,removeactualwashfractions')
    adt=g.input('antibodies','Source anti-biotin/PD-1/TIM3 oligo-antibodies andcellhashingantibodies;wholepolymersunassembled',K.MOLECULAR,state=ReporterComplex(False,False))
    ab_component=g.product('antibody_component','Protein component classes ofsameactualoligo-antibody supplies; mixture membership/stoichiometry unassigned',K.MOLECULAR)
    dna_component=g.product('reporter_dna_component','Copresent DNA tag classes ofsameactualconjugate supplies, includingADT andhash tags; unknownwholepolymers, notoneassertedsingleoligo',K.OLIGO)
    g.react('conjugate_inventory','Component inventory ofactualsourceconjugates, notphysicaldissociation orfreshfreeDNAaddition',S.GENERIC,((adt,R.TEMPLATE),),(ab_component,dna_component),changes_molecule=False)
    tagged=g.product('tagged','Cells bearingglycan/epitope/hashreporters;unboundantibodiesstillpresent',K.MOLECULAR,state=ReporterComplex(True,True))
    g.react('antibody_stain','Source cellhashing/CITE-seq antibodystaining afterlectinwashes',S.BINDING,((x,R.TEMPLATE),(ab_component,R.STRUCTURAL),(dna_component,R.STRUCTURAL)),(tagged,),changes_molecule=False)
    x=clean('antibody_wash',tagged,'Source poststainingcellwashes;counts/prepare~1000cells/uL')
    x=op('gem',x,'Source Chromiumsingle-cell gel-beademulsioncapture/barcoding','Historical5prime VDJFeatureBarcodingkit;exactguide/revisionunreviewed')
    x=op('capture_copy',x,'Source5primeGEX/TCR/featurecapture andcopying continuation;reporterDNAandRNAtemplateclasses distinct','HistoricalRTandcapturecomponents unassembled')
    # Manufacturer-directed operations remain coarse; no 3prime CITE-seq size selection or current kit defaults supplied.
    x=op('source_processing',x,'Manufacturer-directed5prime kitprocessing into source library inputs;exactintermediateoperations unreviewed','HistoricalSingleCellVDJKitwithFeatureBarcoding')
    gene=g.product('gene_input','Source5primeGEXlibraryaliquot;notclaimedpureRNAderivedfraction',K.MOLECULAR)
    tcr=g.product('tcr_input','SourceTCRenrichmentaliquot;maycontainothercDNA',K.MOLECULAR)
    feature_input=g.product('feature_input','Sourcecombined5primeHTO/ADTlibraryinput',K.MOLECULAR)
    g.react('library_allocation','Source GEX,TCR andcombinedADT/HTO preparation interfaces; physicalallocationstage/amounts unassembled',S.MAGIC,((x,R.TEMPLATE),),(gene,tcr,feature_input),reason='Source162 names three actual library preparations. Their allocation/workup schedule is not reconstructed as a literal three-way split or pure-modality fractionation.')
    gene=op('gex_library',gene,'Source5primegeneexpressionlibrarypreparation','Historicalmanufacturerlibraryreagents;unknownmoleculargeometry')
    tcr=op('tcr_enrichment',tcr,'SourcemouseTcellVDJenrichment andTCRlibrarypreparation','ChromiumSingleCellVDJEnrichmentKit,mouseTcell;exacthistoricalrecipeunreviewed')
    f=op('feature_library',feature_input,'Combined5primeHTO/ADTlibrarypreparation andi7indexing','Chromium5primeFeatureBarcodeLibraryKit+i7MultiplexKitNSetA;exactindexoligosunreviewed')
    for x,label in ((gene,'Source5primeGEX'),(tcr,'SourceTCR'),(f,'CombinedsourceHTO/ADT')):
        end(g,x,label+' acquisitionexplicitlyreported;sequencer/readcycles/indexorientation/wholeprimersunassembledinselectedSUGARsource. NoRNAseq75SE orCRISPRscreenNextSeqrunborrowedfromseparatevalidationbranches.')
    return incomplete(g)
def render_page():
    g=protocol_graph();return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_sugar.html">Research notes</a></p>'+''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ('gex_library_product','tcr_enrichment_product','feature_library_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Lectinbinding precedes physicalwash',g.steps['lectin_wash'].inputs[0].material,'lectin_stained')
    c('Antibodystaining useswashedcells',g.steps['antibody_stain'].inputs[0].material,'lectin_wash_kept')
    c('Stainedstate retainsaffinityreporterclass',g.materials['tagged'].state.reporter_bound,True)
    c('CombinedADT/HTOnot inventedseparatelibraries',len(g.steps['library_allocation'].outputs),3)
    c('No RNA-onlyvalidationrunassignment',all(not plan.read_lengths for plan in g.sequencing.values()),True)
    c('ReporterDNA derives from same source conjugate inventory',g.steps['antibody_stain'].inputs[2].material in g.steps['conjugate_inventory'].outputs,True)
    c('Three library names do not assert a physical three-way split',g.steps['library_allocation'].kind,S.MAGIC)
    c('Partialprivategraph',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'sciadv.abe3610_PMC7895430.xml.txt'
    if p.exists():
        t=p.read_text();c('Original explicitcombinedHTO/ADTlibrary','5′ hashtag oligos (HTO)/ADT (one library)' in t,True);c('Original5primeVDJkitratherthan3primeCITEdefault','10× Single Cell V(D)J Kit with Feature Barcoding' in t,True)
    else:c.skip('Originaloperationtranscription','Externaloriginalabsent')
    c.report()
