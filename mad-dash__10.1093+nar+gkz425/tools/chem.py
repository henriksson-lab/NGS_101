"""Published MAD-DASH study route; incomplete molecular geometry stays explicit."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
from chemdraw import Construct,Segment,Scene,feature,revcomp
from protocol_graph import (ProtocolGraph,MaterialKind as K,InputRole as R,
    StepKind as S,SequencingPlan,SourceItem,Terminal,Disposition)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

SOURCE='Published2019 DOI10.1093/nar/gkz425 VORPMC6698666 Methods52–62; publishedSupplementalMaterials_v3 Methods1085–1468 and originalTables_v2 workbookTableS1. Supersededpreprintexcluded.'
A3='ACGGGCTAATATTTATCGGTGG'
A5_CORE=il.TRUSEQ_READ1[-25:]
DIVERSITY=('CAGTCG','TGACTC','GCTAGA','ATCGAT')
RT_PRIMER='GCT'+revcomp(A3)
SAMPLE_INDEX='CAAACGC' # Published seven-base Index01 inventory example; no study sample assigned.
INDEX_FEATURE=feature('mad_dash.sample','sample_index','fixed',note='Published TableS1 Index01 component example, not a sample assignment.')

@dataclass(frozen=True)
class ThirdAdapterChemistry:
    sequence:str=A3
    published_words:str='SI1135–1145 prints adenylation on3primeend and3carbonspacer on5primeend; main says preadenylated without resolving polarity.'
    terminal_assignment:str='unresolved; full authoritative ordering/end chemistry not verified'
    compatible_hypothesis:str='Inferred conventional ligation-compatible5primeadenylate/3primeblock; NOT adopted as a source-established reagent.'

@dataclass(frozen=True)
class RNAAdapterInventory:
    dna_alphabet_body:str
    sugar:str='All source rN residues are RNA; DNA alphabet in coordinate view only'
    terminal_chemistry:str='Full end/modification inventory unassembled; no5primeP assumption'

@dataclass(frozen=True)
class SelectedBreak:
    """Caller-selected strand break; no PAM, cleavage coordinate or efficiency prediction."""
    original:str
    cut:int
    inferred:bool=True
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGT') or type(self.cut)is not int or not 0<self.cut<len(self.original):
            raise ValueError('Explicit DNA and an internal caller-selected cut required')
    @property
    def physical_products(self):return (self.original[:self.cut],self.original[self.cut:])

def body(name,bases):return Construct([Segment(name,bases)])
def indexed_primer():return Construct([Segment('P7',il.P7),Segment('source Index01',SAMPLE_INDEX,feature=INDEX_FEATURE),Segment('source RT arm',RT_PRIMER)])
def rt_binding():
    sc=Scene();sc.strand('adapter_body',[Segment('source A3 ordering body only',A3)],label='Source DNA ordering body; terminal polarity unresolved')
    sc.anneal('rt',[Segment('source RT5prime extension','GCT'),Segment('source binding core',revcomp(A3))],to='adapter_body',pair=('source binding core','source A3 ordering body only'))
    return sc
def missing_scene(label):
    sc=Scene();sc.strand('population',[Segment('unassembled','['+label+']',placeholder=True)]);return sc
def ligation_geometry():
    return AdapterLigationGeometry('unavailable','unavailable',(missing_scene('Source oligo bodies known; complete RNA/substrate/end/junction geometry unassembled'),),'PublishedSI reverses the apparent adenylate/spacer end assignments; no silent repair or full ligation-ready geometry asserted.')

def operation(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete chemistry unassembled'):
    catalytic=kind in {S.DIGESTION,S.DEPROTEINIZATION}
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalytic else K.CONDITION)
    y=g.product(key+'_product',label+'; original molecules and unremoved companions persist',K.MOLECULAR)
    kw={'reason':'Source operation known; full molecular inventory and product geometry unassembled.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalytic else R.CONDITION)),(y,),**kw);return y
def separate(g,key,x,label):
    cond=g.input(key+'_selection',label,K.SELECTION)
    y=g.product(key+'_kept','Actual source-retained fraction: '+label,K.MOLECULAR)
    waste=g.product(key+'_removed','Actual source-separated fraction; imperfect partition/carryover not predicted',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(cond,R.SELECTION)),(y,waste));g.discard(waste,'Source physical separation only, not quantitative selectivity.');return y
def pcr(g,key,x,cycles):
    f=g.input(key+'_f','Published ReversePCRprimer; fresh aliquot of same source species in PCR1/PCR2',K.OLIGO,state=body('source PCR primer',il.NEBNEXT_UNIVERSAL_PRIMER))
    r=g.input(key+'_r','Published ForwardPCRIndex01; same sample index retained in PCR1/PCR2',K.OLIGO,state=indexed_primer())
    enzyme=g.input(key+'_phusion','Source Phusion2xHighFidelity mix',K.CATALYST)
    y=g.product(key+'_product','New PCR descendants plus original templates, cleavage remnants, primers and reagents; no physical deletion of non-amplifying fragments',K.MOLECULAR)
    g.react(key,'Source '+str(cycles)+'cycle '+key,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(missing_scene('Full RNA/RT/PCR/library primer geometry unassembled'),reason='Source primer bodies/cycle counts known; full first-copy/final duplex products unassembled.'))
    return y

def protocol_graph():
    g=ProtocolGraph('MAD-DASH — published brain/multiplex-adapter/bead-cleanup study route',authority=SOURCE,representation='published_instance')
    x=g.input('rna','Source50ng purchased human brain totalRNA; all RNA classes retained, no initial size selection',K.MOLECULAR,primary=True)
    a=g.input('adapter3','Published A3 ordering body; preadenylation known but SI end-polarity wording unresolved',K.ADAPTER,state=body('A3 ordering body only',A3),chemical_state=ThirdAdapterChemistry())
    lig=g.input('ligase2','Source truncated T4RNA ligase2; source main25C1h',K.CATALYST)
    y=g.product('three_ligated','Source3prime-ligation mixture: intended RNAs, largerRNA, free adapter, failed substrates and reagents',K.MOLECULAR)
    g.react('ligate3','Source main3prime-adapter ligation; complete end geometry unresolved',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=ligation_geometry())
    primer=g.input('rt_primer','Source retained RTprimer, added before5prime ligation',K.OLIGO,state=body('RTprimer',RT_PRIMER))
    primed=g.product('primed','Source adapter/RTprimer annealing mixture; includes free/primed adapter and RNA/reagents, no perfect dimer prevention',K.MOLECULAR)
    g.react('anneal_rt','Source main25C1h RTprimer annealing',S.HYBRIDIZATION,((y,R.TEMPLATE),(primer,R.STRUCTURAL)),(primed,))
    adapters=[]
    for i,div in enumerate(DIVERSITY,1):
        seq=A5_CORE+div
        adapters.append(g.input('adapter5_'+str(i),'TableS1RNA5primeadapter'+str(i)+'; terminal diversity is not a barcode/UMI',K.ADAPTER,state=body('Source RNA body, DNA alphabet',seq),chemical_state=RNAAdapterInventory(seq)))
    pool=g.product('adapter5_pool','Source equimolar four-RNA-adapter stock; no preferred member assigned',K.MOLECULAR)
    g.react('pool_adapter5','Source physical four-adapter pool',S.POOL,tuple((a,R.TEMPLATE) for a in adapters),(pool,))
    lig1=g.input('ligase1','Source T4RNAligase1/ATP;adapterstock70C2min thenice before addition',K.CATALYST)
    y=g.product('five_ligated','RNA/adapter ligation mixture includes intended products and adapter dimers, original/failed substrates and reagents',K.MOLECULAR)
    g.react('ligate5','Source main25C1h RNA5prime-adapter ligation',S.LIGATION,((primed,R.TEMPLATE),(pool,R.STRUCTURAL),(lig1,R.CATALYST)),(y,),geometry=ligation_geometry())
    used=g.product('rt_input','Source11ul of mainstudy1:2 split proceeds toRT',K.MOLECULAR)
    rest=g.product('ligation_remainder','Other ligation preparation; mainstudy later fate unassigned',K.MOLECULAR)
    g.react('allocate_rt','Source mainstudy ligation allocation; SI26.5to11ul backup recipe not merged',S.SPLIT,((y,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Mainstudy fate unassigned; alternateSI explicitly stored backup is not silently imposed.'))
    # The original RTprimer is already inside rt_input; do not add another aliquot.
    x=operation(g,'rt',used,'Source SuperScriptII copies from retained RTprimer; originalRNA/freeprimer and unextended products remain',reagent='Source SuperScriptII/RTmix; no new RTprimer aliquot')
    x=pcr(g,'pcr1',x,4)
    x=separate(g,'pre_mad_clean',x,'Source1.8xAMPure cleanup; retain source10ul eluate forMAD-DASH')
    guides=g.input('guide_pool','Source separately IVT-generated four adapter-dimer sgRNAs, normalized/equimolar pooled single-use stocks; whole RNA/scaffold/terminal inventory unassembled',K.MOLECULAR)
    cas=g.input('cas9','Source Cas9 protein binding stock plusNEBuffer3.1/BSA; not SI NEBuffer3/20min variant',K.MOLECULAR)
    rnp=g.product('rnp','Source pre-complexed guide/Cas9 preparation including unbound components',K.CATALYST)
    g.react('precomplex','Source main37C15min Cas9/guide precomplex beforeDNA addition',S.BINDING,((guides,R.TEMPLATE),(cas,R.STRUCTURAL)),(rnp,))
    cut=g.product('cut_mixture','Source conditional Cas9 cleavage products, uncut target/offtarget DNA and bound/unbound RNP; no perfect depletion or physical removal',K.MOLECULAR)
    g.react('cleave','Source main37C2h treatment of limited-PCR dsDNA',S.DIGESTION,((x,R.TEMPLATE),(rnp,R.CATALYST)),(cut,))
    x=operation(g,'rnase',cut,'Source RNaseA37C15min digests susceptible RNA/guide molecules; fragments remain',S.DIGESTION,'Source RNaseA')
    x=operation(g,'protease',x,'Source ProteinaseK37C15min; protein/peptide remnants remain',S.DEPROTEINIZATION,'Source ProteinaseK')
    x=operation(g,'pk_inactivate',x,'Source95C15min inactivatesProteinaseK before immediatePCR2; no prePCR2 cleanup',S.THERMAL_TREATMENT,'Source heating')
    x=pcr(g,'pcr2',x,11)
    x=separate(g,'first_size',x,'Source0.9x45ulbeads; retain transferred supernatant after magnet, high-MW bead fraction removed')
    x=separate(g,'second_size',x,'Source printed90ulbead addition/1.8xoriginalPCRvolume, wash/elute20ul; no reinterpreted effective ratio')
    keep=g.product('production','Source library production retained independently from analytical qPCR aliquot',K.MOLECULAR);qc=g.product('qc','Source KAPA library-quantification analytical aliquot',K.QC)
    g.react('quantify','Source analytical KAPAqPCR allocation',S.SPLIT,((x,R.TEMPLATE),),(keep,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Analytical material does not become sequencing production.'))
    x=operation(g,'normalize',keep,'Source study5nM library normalization; not SI1.25nM recommendation')
    peers=g.input('other_libraries','Other source independently prepared matching-cohort libraries; exact batch identities unassembled',K.MOLECULAR)
    pooled=g.product('sequencing_pool','Source production library pool; exact full strand/read geometry unresolved',K.LIBRARY)
    g.react('pool','Source5nMpool:2ul adapter-dimer-depleted versus1ul undepleted library inputs where included; exact batch membership unassigned',S.POOL,((x,R.TEMPLATE),(peers,R.TEMPLATE)),(pooled,))
    g.sequence(pooled,SequencingPlan(limitation='Actualbrain/bead-cleanup studyHiSeq2500SE50; plasma/bead-cleanup studyNextSeqSE75 separate branch. Seven-base Index01 PCRsampleindex is distinct from six-base RNAadapterdiversity. Exactindex/readprimer/fullend geometry unassembled. Cleaved fragments persist physically; computational filtering is not a molecular discard.'))
    g.add_source_item(SourceItem('pending','Direct Roberts2015/Alon2011 predecessor and historical adapter-order/ligation-end originals','PublishedSI andmain differ in end-polarity wording, reaction formulations/input/normalization and nominal PCR2cleanupvolume; both retained innotes. No full-ready adapter invention, target-efficiency prediction or sourceabsence claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();br=SelectedBreak('ACGTACGT',3)
    c('Selected cleavage conserves all physical bases',''.join(br.physical_products),br.original)
    c('Illustrative selected cut inferred by default',br.inferred,True)
    c.raises('Terminal cuts cannot masquerade as internal breaks',lambda:SelectedBreak('ACGT',0))
    c('RTprimer full binding core derived from source A3',RT_PRIMER[3:],revcomp(A3));rt_binding()
    g=protocol_graph()
    c('PCR1 precedes nuclease treatment',g.steps['cleave'].inputs[0].material,'pre_mad_clean_kept')
    c('Guide digestion does not physically remove cleavage products',g.steps['protease'].inputs[0].material,'rnase_product')
    c('PCR2 follows PKinactivation without invented cleanup',g.steps['pcr2'].inputs[0].material,'pk_inactivate_product')
    c('Remainder is not inferred stored from different SI volume branch',g.terminals['ligation_remainder'].disposition,Disposition.UNRESOLVED)
    c('Analytical QC material excluded from production pool',g.steps['normalize'].inputs[0].material,'production')
    c('Both PCR rounds retain the same sample feature',g.materials['pcr1_r'].state.segments[1].feature,g.materials['pcr2_r'].state.segments[1].feature)
    c('Unresolved source polarity prevents finished publication',bool(g.publication_problems()),True)
    src=Source(sources_root()/'mad-dash__10.1093+nar+gkz425'/'nar_gkz425_supp_gkz425_supplemental_files_Hardigan-et-al_MAD_DASH_Supplemental_Tables_v2.xlsx.txt','Published TableS1 original workbook strings/rows inspected')
    if have(c,src,label='Published TableS1'):
        t=src.path.read_text();c('Source RT and whole PCR primer bodies retained',all(x in t for x in (RT_PRIMER,il.NEBNEXT_UNIVERSAL_PRIMER,il.P7+SAMPLE_INDEX+RT_PRIMER)),True)
        c('Whole four RNA adapter bodies preserve source rN spelling',all(''.join('r'+b for b in (A5_CORE+d).replace('T','U')) in t for d in DIVERSITY),True)
    c.report()

if __name__=='__main__':selftest()
