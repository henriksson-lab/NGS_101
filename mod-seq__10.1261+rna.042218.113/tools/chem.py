"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.'):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),reason=reason);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc


from protocol_graph import Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from illumina import P5,P7
SOURCE='Talkish et al.2014 DOI10.1261/rna.042218.113 publishedPMC3988572 Methods98–116/Fig1; official2019corrigendum DOI10.1261/rna.072405.119 primary publisher indexed text inspected, local XML retrievalHTTP500.'
ORIGINAL_REVERSE_ARM='TGTACTGGAGTTCAGACGTGTGCTCTTCCG'
CORRECTED_REVERSE_ARM='GTGACTGGAGTTCAGACGTGTGCTCTTCCG'
CORRECTION='Official corrigendum2019 DOI10.1261/rna.072405.119 changes originalMethods116 reversePCR arm TGTACTGG… to GTGACTGG…; model uses corrected value. Cached originalPMC main retains old arm; full correction file remains unretrieved, primarypublisher indexed correction text read.'
@dataclass(frozen=True)
class NucleotideBlock:
    name:str
    sequence:str
    def __post_init__(self):
        if not self.name or not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('named exact DNAblock required')
@dataclass(frozen=True)
class CarbonSpacer:
    name:str
    code:str='SpC18'
    def __post_init__(self):
        if not self.name or self.code!='SpC18':raise ValueError('source C18 node required')
@dataclass(frozen=True)
class RNAOriginCopy:
    name:str='cdna'
    origin:str='source RNA-derived DNA extension; selected stop/runoff endpoint unassigned'
@dataclass(frozen=True)
class Backbone:
    nodes:tuple
    circular:bool=False
    def __post_init__(self):
        if not isinstance(self.nodes,tuple) or not self.nodes or not all(isinstance(n,(NucleotideBlock,CarbonSpacer,RNAOriginCopy)) for n in self.nodes) or len({n.name for n in self.nodes})!=len(self.nodes):raise ValueError('ordered unique nucleotide/spacer/copy nodes required')
    def extend(self):
        if self.circular or any(isinstance(n,RNAOriginCopy) for n in self.nodes):raise ValueError('source primer must be linear before RT extension')
        return replace(self,nodes=self.nodes+(RNAOriginCopy(),))
    def close(self):
        if self.circular or not isinstance(self.nodes[-1],RNAOriginCopy):raise ValueError('closure needs extended source linear5phosphateDNA/3OH class')
        return replace(self,circular=True)
    def sequence(self):
        if any(not isinstance(n,NucleotideBlock) for n in self.nodes):raise ValueError('cannot flatten a C18 backbone into DNA sequence')
        return ''.join(n.sequence for n in self.nodes)
    def scene(self):
        sc=Scene();sc.strand('backbone',[Segment(n.name,n.sequence) if isinstance(n,NucleotideBlock) else Segment(n.name,'[C18 non-nucleotide linker]' if isinstance(n,CarbonSpacer) else '['+n.origin+']',placeholder=True) for n in self.nodes],mod5='' if self.circular else 'phosphate',label='source circle opened for display' if self.circular else 'source primer/cDNA');return sc

def rt_primer():return Backbone((NucleotideBlock('five','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC'),CarbonSpacer('spacer1'),NucleotideBlock('middle','CACTCA'),CarbonSpacer('spacer2'),NucleotideBlock('three','TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG')))
@dataclass(frozen=True)
class MixedLinker:
    dna:str='ATCGT'
    rna:str='AGGCACCUGAAA'
    def __post_init__(self):
        if set(self.dna)-set('ACGT') or set(self.rna)-set('ACGU'):raise ValueError('source linker sugar blocks required')
    def scene(self):
        sc=Scene();sc.strand('five_linker',[Segment('DNA',self.dna),Segment('RNA',self.rna,note='source lower-caseaggcaccugaaa is RNA')],label='5prime EndAcceptor17.93R: ATCGTaggcaccugaaa');return sc
    def capture_probe(self):return revcomp((self.dna+self.rna).replace('U','T'))
@dataclass(frozen=True)
class RTMixture:
    dna:Backbone
    def scene(self):
        sc=self.dna.scene();sc.strand('original_rna',[Segment('retained_rna','[modified/end-ligated RNA templates still present before NaOH]',placeholder=True)]);sc.strand('companions',[Segment('carryover','[unextended primers/reagents/remnants]',placeholder=True)]);return sc

def unknown_ligation(statement):
    sc=Scene();sc.strand('unassembled',[Segment('ends','[source population terminal/junction geometry unassembled]',placeholder=True)]);return AdapterLigationGeometry('unavailable','unavailable',(sc,),statement)

def reverse_primer():
    sc=Scene();sc.strand('corrected_reverse',[Segment('p7',P7),Segment('sampleindex','N'*6,placeholder=True,feature=feature('modseq.sample_index','sample_index','unknown')),Segment('correctedarm',CORRECTED_REVERSE_ARM)],label='corrected2019 reversePCR family; study index unassigned');return sc
SOURCE_OLIGO_EXCEPTIONS={'reverse':CORRECTION+' No claim of absence; local main archived pre-correction value.'}

def graph_for_control(control='dms'):
    if control not in {'dms','none','stop'}:raise ValueError('source DMS/noDMS/postquenchDMScontrol required')
    g=ProtocolGraph('Mod-seq — '+control,authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source WT orrpl26-deletion yeast population',K.MOLECULAR,primary=True)
    if control=='dms':x=op(g,'probe',x,'Source100mMDMS30°C2min',S.CHEMICAL_CONVERSION,reagent='freshDMS/ethanol',product='RNA/cell material with chemically methylated RNA classes; stop probabilities/site occupancy not predicted')
    x=op(g,'quench',x,'Source ice and0.6Mbeta-mercaptoethanol/isoamylalcohol quench',S.QUENCHING,reagent='source quenching solution')
    if control=='stop':x=op(g,'postquench_dms',x,'Source STOPcontrol receivesDMS only after quencher addition',S.CHEMICAL_CONVERSION,reagent='sourceDMS afterquencher; validates quenching, no assumed zero residual reaction')
    x=clean(g,'harvest',x,'Source5000g5min pellet/wash, phenolRNAextraction andethanolrecovery')
    x=op(g,'fragment',x,'Source3µgtotalRNA AmbionFragmentationReagents70°C3min',S.FRAGMENTATION,reagent='source historicalfragmentation reagents')
    x=clean(g,'fragment_cleanup',x,'SourceZymoRNAclean>17ntprotocol')
    x=op(g,'dephosphorylate',x,'SourceT4PNK37°C1h phosphatase reaction for fragmentedRNA2prime/3primecyclicphosphates',S.DEPHOSPHORYLATION,catalyst=True,reagent='sourceT4PNK/SUPERaseIn mixture withoutaddedATP; fullenddistribution unassembled')
    x=op(g,'phosphorylate',x,'SourceT4PNK/ATP1mM37°C1h5primephosphorylation',S.PHOSPHORYLATION,catalyst=True,reagent='sourceT4PNK plusATP')
    x=clean(g,'end_cleanup',x,'SourceZymoRNAclean column')
    a=g.input('three_adapter','Source universalmiRNA linker5-rAppCTGTAGGCACCATCAAT-NH2-3; NH2block not an extra nucleotide',K.ADAPTER,state=oligo_scene('three_adapter','CTGTAGGCACCATCAAT',mod5='rApp',mod3='NH2'))
    cat=g.input('rnl2','SourceT4RNAligase2truncated16°Covernight',K.CATALYST);y=g.product('lig3_mix','3prime-adapter-ligated RNA plus unligated RNA/adapter/enzyme remnants',K.MOLECULAR)
    g.react('lig3','Source3prime adapterligation after80°C2min/ice',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=unknown_ligation('Source3prime adenylated17base/NH2adapter; populationend/junction geometry unassembled.'))
    x=clean(g,'lig3_clean',y,'SourceZymo>200ntprotocol')
    a=g.input('five_adapter','Source mixed-sugar5primeEndAcceptor17.93R',K.ADAPTER,state=MixedLinker().scene(),chemical_state=MixedLinker())
    cat=g.input('rnl1','SourceT4RNAligase1/QuickLigasebuffer22°Covernight',K.CATALYST);y=g.product('lig5_mix','5prime-linkered RNA population plus unligated material and source reaction remnants',K.MOLECULAR)
    g.react('lig5','Source5prime linkerligation afterlinker70°C2min/ice',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=unknown_ligation('Source source-known mixedDNA/RNA5linker; complete ligatingends not assembled from predecessor.'))
    x=clean(g,'lig5_clean',y,'SourceZymo>200ntprotocol')
    probe=g.input('rna_capture_probe','Source5-TTTCAGGTGCCTACGAT-3-BiotinTEG',K.OLIGO,state=oligo_scene('rna_probe',MixedLinker().capture_probe(),mod3='Biotin-TEG'))
    bead=g.input('capture_beads','Source streptavidinMyOneC1beads',K.SELECTION)
    bound=g.product('rna_bound','5linker-enriched RNA/bead/probe population, not perfectspecificity',K.MOLECULAR);unbound=g.product('rna_unbound','Separatedcaptureflowthrough; exact composition/fateunassembled',K.MOLECULAR)
    g.react('rna_capture','Sourceprobe0.67µM90°C90s,37°C→RT3°C/min,streptavidincapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(probe,R.STRUCTURAL),(bead,R.SELECTION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Actual removed capturefraction, no downstreamfate inferred.'))
    x=g.product('rna_eluted','Source80°C2minTEmagnet-recovered RNA-supernatant',K.MOLECULAR);beads=g.product('elution_beads','Separatedbead/probe/retainedRNAcompanion, dispositionunassembled',K.MOLECULAR)
    condition=g.input('elution','SourceTE/80°C2min/immediatemagnet',K.CONDITION);g.react('rna_elute','Source thermalRNAelution and magnetseparation',S.CLEANUP,((bound,R.TEMPLATE),(condition,R.CONDITION)),(x,beads));g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'No perfect release or beadfate inferred.'))
    p=g.input('rt_primer','Source5phosphateDNAprimer withtwoC18nodes; specificBrar2012component cited',K.OLIGO,state=rt_primer().scene(),chemical_state=rt_primer());rt=g.input('ssii','SourceSuperScriptII50°C',K.CATALYST)
    mix=RTMixture(rt_primer().extend());y=g.product('rt_mix','Actual RTmixture retaining RNA and spacer-containing cDNA class',K.MOLECULAR,state=mix.scene(),chemical_state=mix)
    g.react('rt','Source reverse transcription',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(p,R.STRUCTURAL),(rt,R.CATALYST)),(y,))
    x=op(g,'hydrolyze',y,'Source100mMNaOH98°C20minRNAhydrolysis',S.CHEMICAL_CONVERSION,reagent='sourceNaOH/heat',product='C18cDNA plusRNAhydrolysisproducts/reagents, not physical removal untilprecipitation/gel')
    x=clean(g,'cdna_size',x,'Sourceethanolprecipitation then8%ureaPAGE150–330bpclass')
    cat=g.input('circligase','SourceCircLigase100U60°C1.5h',K.CATALYST);circle=rt_primer().extend().close();y=g.product('cdna_circle','SelectedcircularcDNAclass retainingC18nodes, plus unclosed/remnantpopulation unquantified',K.MOLECULAR,state=circle.scene(),chemical_state=circle)
    g.react('circularize','SourcecDNAcircularization',S.CIRCULARIZATION,((x,R.TEMPLATE),(cat,R.CATALYST)),(y,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='SourceMethods114; exactpriming/spacertraversalPCRpath unassembled.'))
    probe=g.input('cdna_depletion_probe','SourceBiotinTEG-5-ATCGTAGGCACCTGAAA-3; opposite orientation toRNAcaptureprobe',K.OLIGO,state=oligo_scene('cdna_probe','ATCGTAGGCACCTGAAA',mod5='Biotin-TEG'));bead=g.input('depletion_beads','SourceMyOneC1streptavidinbeads',K.SELECTION)
    kept=g.product('depleted_cdna','Remaining cDNA enriched against5linker-completecopies; no puretrue-stop assertion',K.MOLECULAR);bound=g.product('depletion_bound','Physically removed source5linker-cDNA/probe/beadfraction',K.WASTE)
    g.react('subtractive','Sourceanneal oppositebiotinprobe thenstreptavidin subtractivehybridization',S.AFFINITY_CAPTURE,((y,R.TEMPLATE),(probe,R.STRUCTURAL),(bead,R.SELECTION)),(kept,bound));g.discard(bound,'Source removes5linker-cDNAcontainingfraction from sequencingpreparation; not computationalreadfilter.')
    x=clean(g,'depletion_recover',kept,'Sourceethanolprecipitation ofremainingcDNA')
    f=g.input('forward','SourceP5forwardPCRprimer',K.OLIGO,state=oligo_scene('forward',P5));r=g.input('reverse','Corrected2019 reversePCRfamily; exactsampleindexunassigned. '+CORRECTION,K.OLIGO,state=reverse_primer());poly=g.input('phusion','SourcePhusion/HFbuffer12cycles64°C',K.CATALYST);amp=g.product('amp','PCRproducts plusoriginalspacer-containingtemplates/remnants; source-index/geometryunassembled',K.MOLECULAR)
    g.react('pcr','Source12cyclePhusionPCR usingofficialcorrectedreversearm',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='ExactC18-circlepolymerasetraversal/product/priminggeometryunassembled; no erasedspacers.'))
    x=clean(g,'final_cleanup',amp,'SourceMinElutePCRpurification thenTapeStationQC/pooling')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='Publishedmain116HiSeq2000; acquiredcycles/primergeometryunassembled. Computational5adapterreadfilter is separate from physicalsubtraction. '+CORRECTION))
    g.add_source_item(SourceItem('pending','ActualSupplementalFiles, historicalvendors and specificBrar/Laupredecessors','MainSupplementalMaterial named; EuropePMCzipcontainsfigurealiasesonly. Correctionfile retrievalfailed; exactprimarypublishercorrectiontext inspectedonline. Noabsenceclaim.'))
    return g

def protocol_graph():return graph_for_control()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();link=MixedLinker();circle=rt_primer().extend().close()
    c('mixedsugarRNAacceptor yields complementaryDNAcaptureprobe',link.capture_probe(),'TTTCAGGTGCCTACGAT')
    c('RNAcaptureprobe3biotin is distinct from cDNAdepletion5biotin',graph_for_control().materials['rna_capture_probe'].state.strands['rna_probe'].mod3,'Biotin-TEG')
    c('closure retains two actual non-nucleotideC18nodes',len([n for n in circle.nodes if isinstance(n,CarbonSpacer)]),2)
    try:circle.sequence();bad=False
    except ValueError:bad=True
    c('C18circle cannotbe erasedintoordinaryDNA',bad,True)
    c('postquenchDMScontrol differs fromnoDMS','postquench_dms' in graph_for_control('stop').steps and 'postquench_dms' not in graph_for_control('none').steps,True)
    c('physicalsubtraction removesboundfraction,notremaininglibrary',graph_for_control().steps['subtractive'].outputs,('depleted_cdna','depletion_bound'))
    c('originalandformalcorrectedPCRarms retaineddistinct',ORIGINAL_REVERSE_ARM!=CORRECTED_REVERSE_ARM,True)
    p=Source(sources_root()/'mod-seq__10.1261+rna.042218.113'/'rna.042218.113_PMC3988572.xml.txt','get_sources --slugs mod-seq__10.1261+rna.042218.113;publishedMethods108/112/116')
    if have(c,p,label='sourceoligomodificationguard'):
        text=p.path.read_text();c('source printsNH2block andmixedDNA/RNA5linker','rAppCTGTAGGCACCATCAAT–NH2' in text and 'ATCGTaggcaccugaaa' in text,True)
        c('source prints twoC18backbonenodes','(SpC18)-CACTCA-(SpC18)' in text,True)
    c.report()
if __name__=='__main__':selftest()
