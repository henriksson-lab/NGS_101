"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,additional=()):
    stock=g.input(key+'_stock',reagent,reagent_kind)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,reagent_role))+tuple(additional),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')

@dataclass(frozen=True)
class SelectedLabel:
    dna:str
    sites:frozenset[int]
    stage:str='5hmC'
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or not isinstance(self.sites,frozenset) or any(type(i)is not int or not 0<=i<len(self.dna) or self.dna[i]!='C' for i in self.sites):raise ValueError('explicit DNA and selected 5hmC-bearing cytosines required; not all genomic C')
        if self.stage not in STAGES:raise ValueError('source chemical class required')
    def glycosylate(self):
        if self.stage!='5hmC':raise ValueError('selected native5hmC class needed')
        return replace(self,stage='N3-glucosyl5hmC')
    def click(self):
        if self.stage!='N3-glucosyl5hmC':raise ValueError('azido-glucose required for click')
        return replace(self,stage=CLICKED)

STAGES={'5hmC','N3-glucosyl5hmC','biotin-click-label'}
CLICKED='biotin-click-label'
SOURCE='Associated manuscriptNIHMS295550/PMC3107705 Results79–95/Fig1–4 and publishedSImethods26–137/Fig5/7; DOI10.1038/nbt.1732. ActualonlineMethods/library originals unreviewed/gated.'

def protocol_graph():
    g=ProtocolGraph('hmC-Seal',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source mousecerebellum genomicDNA; genomic5hmC positions/occupancyunassigned',K.MOLECULAR,primary=True)
    x=op(g,'sonicate',x,'Source sonication~100–500bp',S.FRAGMENTATION)
    donor=g.input('udp_azidoglucose','Actual source UDP-6-N3-Glu molecular glycosyl donor; genomic concentration and transferred fraction unassembled',K.MOLECULAR)
    x=op(g,'glyco',x,'Source betaGT transfers6azidoglucose fromUDP6N3Glu',S.CHEMICAL_CONVERSION,reagent='Actual T4 beta-glucosyltransferase catalyst; genomic dose and packaging unassembled',reagent_kind=K.CATALYST,reagent_role=R.CATALYST,additional=((donor,R.STRUCTURAL),),product='SelectedN3-glucosyl5hmCDNA and otherDNA/enzyme/cofactor remnants; no genomicCsubstitution; exactgenomicconditionsunreviewed')
    click_environment=g.input('click_environment','Source copper-free click environment; genomic reaction conditions unassembled',K.CONDITION)
    x=op(g,'click',x,'Source cyclooctynebiotin compound1 copperfreeclick',S.CHEMICAL_CONVERSION,reagent='Actual source cyclooctyne-biotin compound1 molecular reactant; genomic dose and packaging unassembled',reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL,additional=((click_environment,R.CONDITION),),product='Selectedbiotin-N3-5gmCDNA plus unreactedDNA/freecompound1/enzyme/reagent mixture; no disulfide tag inferred from laterSMRTmethod')
    x=op(g,'affinity_recovery',x,'Source biotin-affinityenrichment/recovery of genomicDNA fordeepsequencing',reagent='Actualsource biotin-affinitycapture/recovery original details unreviewed',product='Source recovered5hmC-enrichedDNA class; co-retainedmaterial and separatedcapturefractions remain chemically unassembled. No inventedbeadtype or DTTcleavage')
    x=op(g,'library_boundary',x,'Source deepsequencing librarypreparation; multiple genomiclibraries actualmain',reagent='Actual onlineMethods/genomiclibrary original unreviewed; no currentkitdefaults',product='Source genomicNGSlibrary population; completeadapters/ends/amplification/readgeometry and physicalremnantfates unassembled')
    g.sequence(x,SequencingPlan(limitation='Actual mainFig4 genomicdeepsequencing; exactplatform/cycles/primers unreviewed. Enrichment maps fragment regions, not guaranteed singlebase calls. Main streptavidin-polymerase-stop experiment and HeLa monomericavidin/HRMS validation remain separate.'))
    g.add_source_item(SourceItem('pending','ActualonlineMethods andchemistry/library original','Associatedmain explicitly delegatesonlineMethods; PMC supplementarydata/methodsaliases gated; availablepublisherMOESM6 contains synthesis/figures/data but not inspectedgenomiclibraryrecipe. No absenceclaim.'))
    return g

def selftest():
    from checks import Check
    c=Check();native=SelectedLabel('AACGCC',frozenset({2,5}));tag=native.glycosylate().click()
    c('chemical labeling retains actual DNA bases and selected site identity',(tag.dna,tag.sites),(native.dna,native.sites))
    c('unlabelled cytosine is not automatically marked',4 in tag.sites,False)
    try:native.click();bad=False
    except ValueError:bad=True
    c('click cannot precede azidoglucose installation',bad,True)
    try:SelectedLabel('AACGCC',frozenset({1}));bad=False
    except ValueError:bad=True
    c('modification cannot be installed on declared adenine',bad,True)
    g=protocol_graph()
    c('Glycosyl donor is an actual molecular contributor',g.materials['udp_azidoglucose'].kind,K.MOLECULAR)
    c('betaGT catalyst and UDP donor enter same genomic reaction',tuple((p.material,p.role) for p in g.steps['glyco'].inputs[1:]),(('glyco_stock',R.CATALYST),('udp_azidoglucose',R.STRUCTURAL)))
    c('Click contributor is a molecular reagent rather than enzyme or condition',(g.materials['click_stock'].kind,g.steps['click'].inputs[1].role),(K.MOLECULAR,R.STRUCTURAL))
    c('Actual onlineMethods gate remains uncovered',any("'pending'" in p for p in g.publication_problems()),True)
    c('source-limited graph remains unpublishable',bool(protocol_graph().publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
