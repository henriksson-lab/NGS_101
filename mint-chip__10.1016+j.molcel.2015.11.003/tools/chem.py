from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
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
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='Associated author manuscript PMC4707994 main70–86/134–170, not asserted identical to published final. TableS2 remains gated.'
def adapter_scene():
    sc=Scene();sc.strand('adapter',[Segment('Unassembled source T7-adapter strand','[T7 adapter strand]',placeholder=True)],mod5='C3 spacer')
    return sc

def protocol_graph():
    b=B('Mint-ChIP selected low-MNase pooled chromatin and four assay libraries',AUTHORITY,'Actual native sample cells20ulPBS; titration QC separate/outside production graph');g=b.g
    b.op('digest','20ul2xlysis+MNase/10mMbutyrate; ice20min/37C10min',S.MAGIC,(b.stock('mnase','Source MNase M0247S'),b.cond('lysis','2x100mMTris8/300mMNaCl/2percentTriton/.2percentdeoxycholate/10mMCaCl2')))
    adapter=g.input('adapter_stock','Actual4ul annealed T7 adapter;500uM resuspended/250uM annealed frozenstock;12.5uM source ligation descriptor',K.MOLECULAR,state=adapter_scene(),chemical_state={'sample_index_1':feature('mint.sample','sample_index','unknown',note='Source index#1 sample identity; full bases/segment coordinate unassembled'),'five_prime_block':'Source5primeC3 spacer, non-nucleotide, no guessed other end','components':'T7 promoter/SBS3/sample barcode; duplex whole order unassembled'})
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene(),),'Source double-stranded T7 adapter and 5prime C3 block positive; whole strand order/substrate ends/occupied products unassembled.')
    b.op('ligation','40ul ligationmix/80ul total; RT2h',S.LIGATION,((adapter,R.STRUCTURAL),b.stock('repair_ligase','ActualEndIt1.6ul+FastLink1.6ul enzyme mixes'),b.stock('ligation_nt','Actual1.5mMATP/1xNTPmix source wording; sugar composition unassigned',K.MOLECULAR,R.STRUCTURAL),b.cond('ligation_buffer','1xkitbuffer/30mMEGTA/5mMbutyrate')),geometry=geometry)
    sample=b.op('sample_dilute','Add80ul dilutionbuffer +5mMbutyrate',S.GENERIC,(b.cond('sample_dilution','50mMTris8/150mMNaCl/1percentTriton/50mMEGTA/50mMEDTA/.1percentdeoxycholate'),))
    carrier=g.input('carrier_cells','Actual independent unlabeled carrier cells up to400000 perChIP',K.BIOLOGICAL);b.x=carrier
    b.op('carrier_digest','Carrier lysed/digested concurrently at selected MNase condition',S.MAGIC,(b.stock('carrier_mnase','Actual separately allocated MNase'),b.cond('carrier_lysis','Source same lysis/ice20min/37C10min')))
    b.op('carrier_no_adapter','Carrier receives40ulwater/30mMEGTA/5mMbutyrate, NO T7adapter',S.GENERIC,(b.cond('carrier_stop','Source carrier water/EGTA/butyrate; no ligase/adapters claimed'),))
    carrier=b.op('carrier_dilute','Carrier receives80ul dilutionbuffer',S.GENERIC,(b.cond('carrier_dilution','Source same dilutionbuffer+5mMbutyrate'),))
    other=g.input('other_indexed','Other actually indexed chromatin samples from SAME low-MNase pool; their upstream factories outside selected boundary',K.MOLECULAR)
    pooled=g.product('chromatin_pool','Actual source low-MNase indexed chromatin plus unlabeled carrier pool; high-MNase pool remains separate',K.MOLECULAR)
    g.react('pool_samples','Pool source low-MNase indexed samples then add carrier',S.POOL,((sample,R.TEMPLATE),(other,R.TEMPLATE),(carrier,R.TEMPLATE)),(pooled,))
    names=('h3','h3k4me3','h3k27ac','h3k27me3');parts=tuple(g.product(n+'_chromatin','Actual pooled chromatin assay allocation '+n,K.MOLECULAR) for n in names)
    g.react('split_assays','Split actual common pool for four ChIPs; proportions unassigned',S.SPLIT,((pooled,R.TEMPLATE),),parts)
    libraries=[]
    for name,start,ab in zip(names,parts,('5ulAbcam1791','1ulMillipore07-473','2ulActiveMotif39133','5ulMillipore07-449')):
        b.x=start;b.op(name+'_ab','Source '+ab+' overnight4C',S.BINDING,(b.stock(name+'_abstock','Actual antibody '+ab,K.MOLECULAR,R.STRUCTURAL),b.cond(name+'_abcond','Rotate overnight4C')))
        bound=g.product(name+'_bound','Actual ProteinG retained fraction, occupancy/carryover unassigned',K.MOLECULAR);flow=g.product(name+'_flow','Actual nonretained fraction',K.WASTE)
        g.react(name+'_capture','25ulProteinG Dynabeads2–4h4C',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock(name+'_beads','Actual25ulProteinG',K.SELECTION,R.SELECTION)),(bound,flow));g.discard(flow,'Source magnetic nonretained fraction');b.x=bound
        b.clean(name+'_wash','Source2xRIPA/1x360mMNaClRIPA/1x250mMLiCl/1xTE200ul each')
        b.op(name+'_elute','62C1hDNAelution',S.MAGIC,(b.stock(name+'_pk','Actual.25mgmlProteinaseK'),b.cond(name+'_elution','10mMTris8/1mMEDTA/.1percentSDS/300mMNaCl')))
        b.clean(name+'_dna','1xSPRI >100bp')
        b.op(name+'_ivt','HiScribeT7Quick37C2–16h',S.MAGIC,(b.stock(name+'_t7','ActualT7RNApolymerase mix'),b.stock(name+'_rntp','ActualkitNTPbuffermix',K.MOLECULAR,R.STRUCTURAL)),reason='Positive T7 RNA amplification, full occupied templates/abortive and untagged classes unassembled; unlabeled carrier not magically discarded.')
        b.op(name+'_dnase','1ulkitDNase37C15min after IVT',S.DIGESTION,(b.stock(name+'_dnase_stock','Actual1ulDNase'),))
        b.clean(name+'_rna','Source Silane RNA cleanup')
        b.op(name+'_rt','Source SSIII RT/50uM RTprimer randomhexamer+SBS12 region',S.MAGIC,(b.stock(name+'_rtmix','Actual SSIII FirstStrandSynthesisSuperMix; enzyme/nucleotide composition unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock(name+'_rtprimer','Actual50uM RTprimer; randomhexamer priming, NOT UMI; whole order unassembled',K.MOLECULAR,R.STRUCTURAL)))
        b.clean(name+'_cdna','1xSPRI >100bp')
        for number in (1,2):
            if number==2:b.split(name+'_pcr2_input','Actual2–40ng input from PCR1 DNA stock','Actual2–40ng PCR2 input','UntransferredPCR1DNA stock; fate unassigned')
            b.op(name+'_pcr'+str(number),'Source PCR'+str(number)+' up to12cycles; usually '+('6'if number==1 else'8'),S.MAGIC,(b.stock(name+'_pcrmix'+str(number),'ActualPfuUltraIIFusionHS master mix aliquot; full catalyst/nucleotide composition unassembled',K.MOLECULAR,R.STRUCTURAL),b.stock(name+'_primers'+str(number),'Actual.2uM forward/reverse'+str(number)+('; reverse2 assay index#2 distinguishes antibody, NOT UMI'if number==2 else'')+'; whole orders unassembled',K.MOLECULAR,R.STRUCTURAL,chemical_state={'assay_index_2':feature('mint.assay.'+name,'feature_barcode','unknown',note='Source reverse2 assay index#2 identifies '+name+' assay; sequence/coordinate unassembled')}if number==2 else{})))
            b.clean(name+'_postpcr'+str(number),'1xSPRI >100bp'if number==1 else'DoubleSPRI.6x/.95x200–700bp')
        b.split(name+'_qc','Source Qubit/BioAnalyzer/KAPAquant analytical aliquots','Actual production library after QC allocation','Actual consumed/retained analytical QC aliquots; quantities unassigned',disposition=Disposition.QC)
        b.split(name+'_pool_input','Select actual assay-library pool contribution, source quantities unassigned','Actual aliquot contributing to final assay pool','Actual untransferred postQC assay-library stock; fate unassigned');libraries.append(b.x)
    final=g.product('assay_library_pool','Actual indexed assay-library sequencing pool; pooling amounts unassembled',K.LIBRARY)
    g.react('pool_assay_libraries','Pool assay-indexed libraries for paired-end sequencing',S.POOL,tuple((x,R.TEMPLATE)for x in libraries),(final,))
    g.sequence(final,SequencingPlan(limitation='NextSeq500/HiSeq2500 paired-end; sampleindex#1 versus assayindex#2; exact runmask/whole oligos/landing unassembled. Randomhexamer is priming class, not dedup UMI.'))
    g.add_source_item(SourceItem('identifier_scope','Source main72–86/146–170 index#1 sample identity before pooledChIP versus index#2 antibody assay at PCR2.','Full TableS2 and exact sequence/feature coordinates unreviewed, not confirmed absent.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Carrier factory is independent',g.steps['carrier_digest'].inputs[0].material,'carrier_cells');c('Sample pool follows source ligation',g.steps['pool_samples'].inputs[0].material,'sample_dilute_product');c('Four actual assay allocations',len(g.steps['split_assays'].outputs),4);c('Final actual assay pool',len(g.steps['pool_assay_libraries'].inputs),4);c('DNase after IVT',g.steps['h3_dnase'].inputs[0].material,'h3_ivt_product');c('PCR2 consumes selected stock allocation',g.steps['h3_pcr2'].inputs[0].material,'h3_pcr2_input_selected');c('Pool uses selected assay aliquot',g.steps['pool_assay_libraries'].inputs[0].material,'h3_pool_input_selected');c('C3 block source-known',g.materials['adapter_stock'].chemical_state['five_prime_block'].startswith('Source5primeC3'),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2015.11.003_PMC4707994.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2015.11.003')):
        t=p.read_text();c('Original carrier receives water','instead of ligation mix' in t,True);c('Original two separate MNase pools','pool 1' in t and 'pool 2' in t,True)
    c('Privatepartial unready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
