"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry,reported_adapter_end

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
 st=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.OLIGO if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; all substrates/products/companion species and reaction remnants retained until stated workup',K.MOLECULAR)
 kw={'reason':'Known source operation; exact molecular states, ends and geometry not yet assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(st,role)),(y,),**kw);return y

def cleanup(g,key,x,label):
 st=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; recovery/carryover not presumed perfect',K.MOLECULAR);w=g.product(key+'_waste','Actual separated physical workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(st,R.CONDITION)),(y,w));g.discard(w,'Source physical separation, not computational filtering.');return y

def split(g,key,x,labels):
 ys=[g.product(key+'_'+str(i),lab,K.MOLECULAR) for i,lab in enumerate(labels)];g.react(key,key,S.SPLIT,((x,R.TEMPLATE),),ys);return ys

def assay(g,key,x,label):
 y,a=split(g,key,x,('Retained preparation after source assay allocation; assay amount unassigned',label));g.mark_terminal(Terminal(a.id,Disposition.QC,'Separate analytical material; no implicit return to preparation.'));return y

def pcr(g,key,x,label,primers):
 ps=[g.input(key+'_p'+str(i),lab,K.OLIGO,state=state) for i,(lab,state) in enumerate(primers)];e=g.input(key+'_enzyme',label+' source catalyst',K.CATALYST);y=g.product(key+'_product','New DNA copies plus original DNA/free primers/reagents',K.MOLECULAR);geom=PCRGeometry.unavailable(unknown('Whole operational PCR geometry not assembled'),reason='Source primer identities do not by themselves establish all inward sites/products.',symmetric_species=len(ps)==1);g.react(key,label,S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in ps),(e,R.CATALYST)),(y,),geometry=geom);return y

def oligo_scene(name,seq,mod5='',mod3=''):
 s=Scene();s.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return s

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

from dataclasses import dataclass
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
SOURCE='Version-qualified PMC6697087.1 author manuscript, overview90–94 and Methods210–248. Selected fourplex RNA-transcript assay; actual Tables S1/S2 are gated, not claimed absent.'
@dataclass(frozen=True)
class UMIComposition:
 four_way_positions:int=20
 two_way_positions:int=9
 @property
 def length(self):return self.four_way_positions+self.two_way_positions
 @property
 def capacity(self):return 4**self.four_way_positions*2**self.two_way_positions
# Composition is source-reported; no synthetic ordering or full primer is supplied.
UMI_COMPOSITION=UMIComposition()
def event_identifier():
 s=Scene();s.strand('isolated_event_region',[Segment('UEI only; full neighboring sites unassembled','N'*20,placeholder=True,feature=feature('dna_microscopy_uei','umi','random'))]);return s

def protocol_graph():
 g=ProtocolGraph('DNA microscopy — selected fourplex transcript preparation',authority=SOURCE,representation='published_instance')
 x=g.input('cells','Source cultured fixed-position cells containing RNA transcripts; selected ACTB beacon and GFP/RFP/GAPDH targets, not genomic-DNA assay',K.BIOLOGICAL,primary=True)
 x=op(g,'fix_permeabilize',x,'Source 4% formaldehyde PBS 15min RT; PBS washes, 0.25% Triton 10min, PBS, 0.1N HCl 2–3min; stock kept 4C')
 x=cleanup(g,'pre_rt_washes',x,'Source physical PBS washes and ddH2O rinse/aspiration before RT')
 rt=g.input('rt_ultramers','Source UMI-bearing RT ultramers, Table S1 fourplex; 850nM each samples1–2 or100nM others; exact sites and 20N+9 two-way position ordering unassembled',K.OLIGO)
 enzyme=g.input('rt_reagents','Source SSIII10U/ul, 400uM dNTP, SuperaseIn1U/ul, 1x buffer, printed4uM DTT',K.CATALYST)
 y=g.product('rt_product','Source beacon/target UMI-tagged cDNA plus original RNA/cells/free RT oligos/reagents',K.MOLECULAR);g.react('reverse_transcription','Source mixture incubated60C3min then42C1h then4C; no reordered preheat',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(rt,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,));x=cleanup(g,'post_rt_washes',y,'Source aspiration and three PBS washes, overnight4C in final wash, ddH2O rinse/aspiration')
 x=op(g,'exoi',x,'Source40ul ExoI1.4U/ul in1x buffer,37C40min',S.DIGESTION,'Source exonuclease I; detailed ssDNA accessibility/retention not predicted')
 x=cleanup(g,'post_exoi_wash',x,'Source three physical PBS washes after ExoI')
 x=op(g,'hydrogel',x,'Source separately prepared PEG components mixed immediately,20ul added,10.8% hydrogel polymerizes1h RT; 4arm acrylate64ug/ul and2arm thiol44ug/ul')
 x=pcr(g,'in_situ_overlap',x,'Source Platinum Taq hydrogel fourplex:95C2min;10x95C30s/68C1min;2x95C30s/55C30s/68C1min;16x95C30s/60C30s/68C1min;68C1min. Initial10 intended linear; UEI overlap-extension, not predicted yield/diffusion',tuple((name,unknown(name+' exact TableS2 sequence unassembled')) for name in ('OE1a','OE4b','psbs12s','s8B','LF primer collection')))
 x=op(g,'dissolve',x,'Source stored−20C, thaw ice;4ul460mM KOH/100mM EDTA/40mM DTT2h ice;72C5min;mix10times')
 x=op(g,'neutralize',x,'Source4ul neutralization stock made400ul1N HCl plus600ul1M TrisHClpH7.5; immediate mixing')
 x=op(g,'protease',x,'Source11.1ul0.35%Tween20/0.35mgml proteinaseK/10mMTrispH8;50C25min',S.DEPROTEINIZATION,'Source proteinase K; no RNase or DNA-specific cleanup invented')
 x=op(g,'dilute',x,'Source55ul10mM TrisHClpH8 added and mixed')
 used,rest=split(g,'plate_transfer',x,('Source85ul transferred to new PCR plate','Untransferred post-dissolution sample; amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'85ul is not assumed to consume whole nominal preparation.'))
 x=cleanup(g,'ampure1',used,'Source0.65x AMPure XP10min;two70%EtOH washes;elute35ul10mMTrispH8')
 used,rest=split(g,'ex_situ_allocation',x,('Source product aliquot diluted1:2 into PCR; exact volume unassigned','Unassigned unused35ul eluate stock if any'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No invented whole35ul consumption or exact remaining volume.'))
 x=pcr(g,'interference_pcr',used,'Source Platinum Taq HiFi0.02U/ul,1.5mM MgSO4,200uM dNTP;95C2min;20x95C30s/68C2min; interference primers intended to suppress new concatemers',tuple((name,unknown(name+' TableS2 sequence unassembled')) for name in ('SBS3LC300nM','rev-ill-214300nM','10T-0E-P3.3uM','10T-OEc-P3.3uM')))
 x=cleanup(g,'ampure2',x,'Source0.65x AMPure XP workup;elute40ul10mMTrispH8')
 used,rest=split(g,'barcode_allocation',x,('Source10ul40ul eluate to barcoding PCR','Source untransferred30ul eluate stock; fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Physical untransferred stock; no computational discard.'))
 x=pcr(g,'barcode_pcr',used,'Source forill-sbs3/rev-ill-X300nM each,Platinum Taq HiFi0.02U/ul,2mM MgSO4,200uM dNTP;95C2min;5x95C30s/58C30s/68C2min;1–5x95C30s/68C2min; exact sample-cycle count unassigned',(('forill-sbs3',unknown('TableS2 primer')),('rev-ill-X carrying fixed8nt SAMPLE index',unknown('TableS2 full primer/index value'))))
 x=cleanup(g,'ampure_final',x,'Source final0.7x AMPure XP purification as above; final elution volume not assigned here')
 return finish(g,x,'Source NextSeq550 PE: SBS3 canonicalR1 103bp and SBS12 canonicalR2 57or151bp, sample-dependent TablesS1/S5. Fixed8nt sample index 5prime-adjacent SBS12 is separate from beacon/target UMIs and20nt UEI. Exact molecular positions/primer landing geometry unassembled. UMI/UEI clustering and spatial image inference are computational, not physical DNA removal. 24plex alternative has distinct primers/program, not merged.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('All three aliquot boundaries retain physical remainder',all(k in g.terminals for k in ('plate_transfer_1','ex_situ_allocation_1','barcode_allocation_1')),True);c('No invented whole-library/read geometry',bool(g.publication_problems()),True);c('Event region has molecular-event meaning',event_identifier().strands['isolated_event_region'].segs[0].feature.role,'umi');c('UMI composition is not 29 unconstrained Ns',UMI_COMPOSITION.capacity<4**UMI_COMPOSITION.length,True)
 src=Source(sources_root()/'dna-microscopy__10.1016+j.cell.2019.05.019'/'j.cell.2019.05.019_PMC6697087.xml.txt','AMmain overview90–92 and Methods230/248')
 if have(c,src,label='SourceAM'):
  t=''.join(src.path.read_text().split());c('Whole canonical sequencing primer bodies verified',all(s in t for s in (TRUSEQ_READ1,TRUSEQ_READ2)),True);c('UEI20nt explicit in overview, not inferred diversity','new20-nucleotideUEI' in t,True);c('Mixed-degeneracy UMI composition verified','ℓ4=20andℓ2=9' in t,True)
 c.report()
if __name__=='__main__':selftest()
