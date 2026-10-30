"""Source-local uvCLAP preparation and tagged-window declarations; no complete library."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Segment,Scene,feature
import illumina as il
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
SOURCE='PublishedPMC5861125 Results45–53/Methods171/183–201; publishedSI/MOESM1 FigS1originalp2 andTable7originalp12. Main/SI two originals; no predecessor recipe imported.'
TAG='AGTCA'
P5_PRIMER=il.TRUSEQ_P5_FULL
P7_PATTERN=il.P7+'N'*6+il.TRUSEQ_READ2

@dataclass(frozen=True)
class TagLayout:
    """Source pattern membership; no observed library/specimen/complete oligo assignment."""
    five_tag:str=TAG
    replicate_pattern:str='DRYYR'
    inferred:bool=True
    def __post_init__(self):
        if len(self.five_tag)!=5 or set(self.five_tag)-set('ACGT') or self.replicate_pattern not in {'DRYYR','DYRRY'}:raise ValueError('Explicit5baseconditiontag and one of two source replicate-pattern families required')
    @property
    def replicate_meanings(self):
        # Both meanings are structured locally; one Segment feature alone cannot express both.
        return (feature('uv_replicate','sample_index','degenerate'),feature('uv_rep_random','umi','degenerate',note='Replicate family membership and within-family variation both carry information; complete read geometry unassembled'))
    def five_window(self):
        return Construct([Segment('UMI left','NNN',placeholder=True,inferred=self.inferred,feature=feature('uv_umi_left','umi','random',group='uv_umi',part='left')),Segment('condition tag',self.five_tag,inferred=self.inferred,feature=feature('uv_condition','sample_index','fixed')),Segment('UMI right','NN',placeholder=True,inferred=self.inferred,feature=feature('uv_umi_right','umi','random',group='uv_umi',part='right'))])

def unknown(label):
    s=Scene();s.strand('mixture',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,reagent='Source catalysts/conditions; complete mixture inventory unassembled',kind=S.MAGIC):
    is_enzyme=kind in {S.DIGESTION,S.DEPROTEINIZATION,S.REVERSE_TRANSCRIPTION}
    stock=g.input(key+'_stock',reagent,K.CATALYST if is_enzyme else K.CONDITION);y=g.product(key+'_product',label+'; companions/reagents retained until source separation',K.MOLECULAR)
    kwargs={'reason':'Source-known operation; exact molecular geometry and full conditions remain unassembled.'}if kind==S.MAGIC else {}
    inp=[(x,R.TEMPLATE),(stock,R.CATALYST if is_enzyme else R.CONDITION)]
    if kind==S.REVERSE_TRANSCRIPTION:inp.append((g.input(key+'_primer','Source barcodedRTprimer; wholeoligo/modifiedendsunassembled',K.OLIGO),R.STRUCTURAL))
    g.react(key,label,kind,inp,(y,),**kwargs);return y

def separate(g,key,x,label):
    r=g.input(key+'_selection',label,K.SELECTION);y=g.product(key+'_retained','Source retained physical preparation; exact composition/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Actual source separated/unbound/matrix fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(r,R.SELECTION)),(y,w));g.discard(w,'Source physical separation only; no perfect target/background classifier.');return y

def protocol_graph():
    g=ProtocolGraph('uvCLAP — source His/streptavidin path with jointly amplified libraries',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source taggedRBP-expressing cells;tag installation/expression upstreamunassembled',K.BIOLOGICAL,primary=True)
    x=op(g,'uv',x,'Main UV-C~254nm printed0.15mJ/cm2; no silent unit conversion to0.15J/cm2')
    x=op(g,'harvest_lysis',x,'Sourcepellet/snapfreeze/thaw/0.5mlHislysis/mildsonication; exactmixture inventory unassembled')
    x=separate(g,'his_capture',x,'SourceHisTagDynabeads10min/lysisbufferwash;FLAG1h branch unselected')
    x=separate(g,'imidazole_release',x,'Source250mMimidazoleeluate retained; originalHisbeadmatrix separated')
    x=separate(g,'streptavidin',x,'SourceMyONEC1 secondcapture/highstringencySDS/NaCl/LiDS/LiCl washes;proteinbiotin is not DNAUMI')
    x=op(g,'rna_trim',x,'SourceNDB/TURBODNase/RNaseI37C3min→ice',reagent='Source2ulTURBODNase plus10uldilutedRNaseI1:2000–1:8000from100Uul in1mlNDB',kind=S.DIGESTION)
    x=op(g,'pnk',x,'LowpHPNK cyclicphosphate-enddephosphorylation37C20min; no known perfect conversion assumed')
    a=g.input('adapter3','Sourcepreadenylated3prime RNA linker perSI FigS1; fullbody/endsunassembled',K.ADAPTER);e=g.input('rna_ligase','MainT4RNAligase1; sourceconditions/unreviewedfullgeometry',K.CATALYST);y=g.product('rna_adapter_mixture','Adapterligated/unligated RNA-protein/reagentmixture',K.MOLECULAR)
    g.react('rna_adapter','SourceRNA3primeadapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown('Full3primeRNA/linker endgeometry unassembled'),),'Pre-adenylated linker source-known,complete end/junction geometry not assembled.'))
    x=separate(g,'adapter_wash',y,'Source wash removes excess adapters')
    x=op(g,'pk',x,'SourceProteinaseK release of3taggedRNA',reagent='SourceProteinaseK/proteolysisconditionsunassembled',kind=S.DEPROTEINIZATION)
    x=separate(g,'column',x,'SourceZymoDNA CleanandConcentratorcolumn recovery')
    x=op(g,'rt',x,'SourceSuperScriptIII/barcodedRT, mixture retained beforepool',reagent='SourceSSIII/dNTPs;exactbuffer/timeunassembled',kind=S.REVERSE_TRANSCRIPTION)
    other=g.input('other_cdnas','Independent biologicalreplicate/condition/mockcontrolcDNAs; individualpreparations unassembled',K.MOLECULAR,primary=True);pool=g.product('cdna_pool','Source pooledcDNAs beforegel andjointPCR; no independentnormalization asserted',K.MOLECULAR);g.react('joint_pool','MixrelevantRTproducts before gel fractionation',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
    x=separate(g,'urea_page',pool,'Source6percent/6MureaPAA fractionatescDNA; exact windows/fractionmembershipunassigned')
    x=op(g,'circle',x,'SourceCircLigasecDNAcircularization;fullmolecules/junction/reagentinventoryunassembled')
    x=op(g,'linearize',x,'Source restrictionlinearization; sourceenzyme/site/cutgeometry not assigned',reagent='Source restrictionenzymestock;identityunassignedinthispass',kind=S.DIGESTION)
    f=g.input('pcr_forward','PublishedScriptSeqforwardprimer',K.OLIGO,state=Construct([Segment('sourceScriptSeqP5',P5_PRIMER)]));r=g.input('pcr_reverse','PublishedScriptSeqindexedprimerpattern; sixfixedsampleindexbasesunassigned,notUMI',K.OLIGO,state=Construct([Segment('sourceP7',il.P7),Segment('PCRsampleindex','N'*6,placeholder=True,feature=feature('uv_pcr_sample','sample_index','fixed')),Segment('sourceRead2handle',il.TRUSEQ_READ2)]));e=g.input('pcr_catalyst','SourcePCRcatalyst/cycles notassignedinthispass',K.CATALYST);lib=g.product('library','Sourcejointlyamplifiedsequencingpreparation, originalDNA/primers/reagents remainunassembled',K.LIBRARY)
    g.react('joint_pcr','SourcejointPCR of sizefractionated/circularized/linearized cDNA',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(e,R.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(unknown('FullsourcecDNA/RT/circle/linearization/PCRgeometryunassembled'),reason='CanonicalPCRbodies/patternknown; no wholelibraryconstructor imported.'))
    g.sequence(lib,SequencingPlan(limitation='Sourcepairedend default; exactcycles/instrument/run/indexassignmentunassembled. Condition/replicate/sizefractiontags distinct. UMIeventdedup/backgroundpeakfiltering computational, not a preparation reaction.'))
    g.add_source_item(SourceItem('pending','FulloperationalRT/linker/circle/cut/PCR/run and sourcecoverage','Publishedselectedscope only; otherfiles/experiments/linkedpredecessors/vendororiginalsunread,notconfirmedabsent.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();layout=TagLayout();parts=layout.five_window().segments
    c('SplitUMIparts sharelogicalgroup',parts[0].feature.group==parts[2].feature.group,True);c('SplitUMIparts retain distinctidentities',parts[0].feature.id!=parts[2].feature.id,True);c('Prototypewindowinferredbydefault',all(s.inferred for s in parts),True);c('Dualreplicatemeanings retained',set(f.role for f in layout.replicate_meanings),{'sample_index','umi'});c('AllRTclassespoolbeforefractionation',g.steps['urea_page'].inputs[0].material,'cdna_pool');c('No fullread/chemicalvalidationclaim',bool(g.publication_problems()),True)
    try:TagLayout(replicate_pattern='NNNNN');bad=False
    except ValueError:bad=True
    c('Unsupportedreplicatepatternrejected',bad,True)
    src=Source(sources_root()/'uvclap__10.1038+s41467-018-03575-4'/'s41467-018-03575-4_PMC5861125.xml.txt','get_sources --slugs uvclap__10.1038+s41467-018-03575-4;publishedMethods187')
    if have(c,src,label='PublishedMain'):
        text=''.join(src.path.read_text().split());c('BothwholepublishedPCRpatterns',P5_PRIMER in text and P7_PATTERN in text,True);c('LiteralUVunit preserved','0.15mJ/cm2'in text,True)
    c.report()
if __name__=='__main__':selftest()
