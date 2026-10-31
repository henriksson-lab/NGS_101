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

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,extra_inputs=()):
    stock=g.input(key+'_stock',reagent,K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CONDITION))+tuple(extra_inputs),(out,),**kw);return out

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

    def reduce(self):
        if self.stage!='biotin-S-S-label':raise ValueError('DTT release requires source disulfide tag')
        return replace(self,stage='HS-label-remnant')
STAGES={'5hmC','N3-glucosyl5hmC','biotin-S-S-label','HS-label-remnant'}
CLICKED='biotin-S-S-label'
SOURCE='Associated manuscriptNIHMS465877/PMC3646335 Methods68–84/94 and publishedSupplementMOESM18 Fig1/Table3; DOI10.1038/nmeth.1779. HistoricalSMRTbellpredecessor unreviewed.'

def protocol_graph():
    g=ProtocolGraph('5hmC chemical-label SMRT sequencing',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourcemESC E14Tg2A genomicDNA;5hmCsites/occupancy not instantiated',K.MOLECULAR,primary=True)
    x=op(g,'sonicate',x,'Source sonication100–500bp',S.FRAGMENTATION)
    gt=g.input('label_gt','Actual genomic-labeling betaGT2uM contribution',K.CATALYST)
    sugar=g.input('label_sugar','Actual genomic-labeling UDP6N3Glc100uM contribution',K.MOLECULAR)
    x=op(g,'label',x,'Source200µgDNA/UDP6N3Glc100µM/betaGT2µM37°C1h',S.CHEMICAL_CONVERSION,reagent='Source50mMHEPESpH7.9/25mMMgCl2',product='SelectedN3-glucosyl5hmCDNA plus unreactedDNA/enzyme/cofactor remnants',extra_inputs=((gt,R.CATALYST),(sugar,R.STRUCTURAL)))
    x=cleanup(g,'label_clean',x,'SourceBioRadMicroBioSpin6 desalting')
    linker=g.input('click_linker','Actual source DBCOdisulfidebiotin150uM contribution',K.MOLECULAR)
    x=op(g,'click',x,'Source150µM DBCOdisulfidebiotin37°C2h',S.CHEMICAL_CONVERSION,reagent='Source37C2h incubation after actual click-reagent addition',product='Sourcebiotin-S-S-N3-5gmC DNA classes plus freeclickreagent/unreactedDNA; no baseconversion',extra_inputs=((linker,R.STRUCTURAL),))
    x=cleanup(g,'click_clean',x,'SourceQiagennucleotideremoval purification')
    beads=g.input('beads','Source100µlMyOneStreptavidinC1beads/B&W0.01%Tween20',K.SELECTION);bound=g.product('bound','Sourcebead-associated labeledDNA plus unquantified nonspecific companions',K.MOLECULAR);flow=g.product('flowthrough','Source separated captureflowthrough/washes',K.WASTE)
    g.react('capture','SourcebiotinDNA binding and washing',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.discard(flow,'Actual removed binding/wash fraction; no perfect specificity asserted.')
    reductant=g.input('dtt_reductant','Actual source DTT50mM in100ul release solution',K.MOLECULAR)
    y=op(g,'dtt',bound,'Source50mMDTT100µl25°C2hgentlerotation',S.CHEMICAL_CONVERSION,reagent='Source25C2hgentlerotation release programme',product='ReleasedHS-N3-5gmCDNA plus beads with residualbiotin/linker products and DTT/reactionremnants; no physical separation yet',extra_inputs=((reductant,R.STRUCTURAL),))
    stock=g.input('magnet','Source magneticsegregation',K.CONDITION);x=g.product('released','SourcecollectedDNA solution retainsHS-labeladduct; no restoration to native5hmC',K.MOLECULAR);beadrest=g.product('retained_beads','Source magnet-separatedbeads/biotin-linkerremnants and possible retainedDNA',K.WASTE)
    g.react('magnet_collect','Sourcecollect releasedDNA supernatant',S.CLEANUP,((y,R.TEMPLATE),(stock,R.CONDITION)),(x,beadrest));g.discard(beadrest,'Actual magnetic bead separation; recovery incomplete/unspecified.')
    x=cleanup(g,'desalt',x,'SourceBioRadMicroBioSpin6 removesDTT')
    x=cleanup(g,'minelute',x,'SourceQiagenMinElutePCRpurification/EBelution')
    x=op(g,'smrtbell',x,'Sourcegenomic A-tailing SMRTbellpreparation perTravers2010ref14',reagent='Source historicalA-tailingSMRTbellmethod; originalrecipe/endgeometryunreviewed',product='Source topologically circular genomicSMRTbell population retainingmodifiedinsert; exacthairpin/junctions/remnants unassembled')
    g.sequence(x,SequencingPlan(limitation='SourcePacBioRS smallinsert acquisition; native retained adduct modulates IPD, not a C-to-T conversion. Historical polymerase/primer geometry unreviewed. Syntheticfourbaseoverhang/exoIII/VII branch is not assigned to genomic library. Computational IPD thresholds/CCS do not alter physical DNA.'))
    g.add_source_item(SourceItem('pending','HistoricalgenomicSMRTbellA-tailing original and PacBioRSstandardprotocol','Main84 explicitly citesTravers2010 DOI10.1093/nar/gkq543; not followed in cappedpass. No ordinary modernkit recipe imported.'))
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
    c('source-limited graph remains unpublishable',bool(protocol_graph().publication_problems()),True)
    released=tag.reduce();c('DTT leaves a DNA-bound thiol-linker remnant',released.stage,'HS-label-remnant')
    c('DTT leaves all nucleotide positions intact',released.dna,native.dna)
    try:native.reduce();bad=False
    except ValueError:bad=True
    c('no reduction-release from native or non-disulfide class',bad,True)
    g=protocol_graph()
    c('Labeling consumes actual source betaGT catalyst',g.steps['label'].inputs[2].material,'label_gt')
    c('DTT release consumes physical reducing reagent before magnetic separation',g.steps['dtt'].inputs[2].material,'dtt_reductant')
    c('Genomic historical SMRTbell scope remains uncovered','pending' in g.source_coverage,False)
    c('beads are physically removed after chemicalreduction',g.steps['magnet_collect'].inputs[0].material,'dtt_product')
    c.report()
if __name__=='__main__':selftest()
