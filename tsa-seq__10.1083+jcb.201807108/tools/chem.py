"""Chen2018 TSA-seq: native pulldown and matched-input chemistry."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, feature
from protocol_model import role
from endprep import EndPreparedDuplex, dA_tailed_scene
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp
TITLE='TSA-seq — distance mapping from a nuclear compartment'
NOTES='01_tsa_seq.html'
SOURCE='<a href="https://doi.org/10.1083/jcb.201807108">Chen et al., Journal of Cell Biology2018</a>, published Methods222–260.'
SUMMARY='Antibody-localized HRP produces a tyramide-biotin gradient in fixed nuclei. Sonicated DNA is split into10%matched input and90%streptavidin pulldown; both enter source KAPA HyperPrep preparation.'
CAVEAT='Distance is inferred from population enrichment, not encoded in a covalent barcode or DNA junction. Physical6nt adapter barcodes do not establish index position/count or6cycle acquisition. The source reports101sequencing cycles and100nt output reads without assigning SE/PE mode here. Native adapter strands and PCR geometry remain unassigned. Condition1labelling and latermodifiedall-at-once capture are selected source procedural references; sample-specific capturevariant pairing remains unassigned.'
SAMPLE_INDEX=feature('tsa_sample_index','sample_index','unknown',note='Source indexed adapters have6nt barcodes; native position/count/member unassigned')
INDEX_ROLE=role('source sample barcode; native position unassigned',length_bp=6,molecular_feature=SAMPLE_INDEX)
FINAL_LIBRARY=Construct([role('native KAPA-prepared endA'),role('genomic DNA copy; source biotin position not a sequencing barcode'),INDEX_ROLE,role('native opposite/index-associated end; architecture unassigned')],name='TSA native indexed library reference')
SEQ_PRIMERS=()
FINAL_CAPTION='Native library roles; sample barcode6nt is source-positive, complete adapter architecture remains unassigned.'
REPORTED_CYCLES=101
REPORTED_READ_NT=100
SEQUENCING_INTRO=f'The paper reports{REPORTED_CYCLES}cycles on HiSeq2000/2500/4000 and{REPORTED_READ_NT}nt reads. SE/PE mode and native index acquisition remain unassigned; pulldown/input enrichment supplies the distance signal.'
SEQUENCING_UNAVAILABLE='Native sequencing-primer sites, read mode and index acquisition remain unassigned.'
SEQUENCING_ENDING=SEQUENCING_UNAVAILABLE
SHOW_READ_CONTENT=False
READ_SPECS=tuple(sp.read_spec(r,'unassigned','Source101sequencing cycles/100nt reads; role/mode unassigned' if r.startswith('Read') else 'Physical6nt sample barcode does not assign index acquisition','undetermined') for r in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
READ_SPECS_NOTE=SEQUENCING_INTRO

def tsa_library():return FINAL_LIBRARY,()
def cellular(stage):
 sc=Scene()
 sc.strand('proteins',[role('nuclear proteins/antibody-localized HRP: '+stage)],label='protein assemblies; independent of DNA strand geometry')
 for label in ('nearer chromosomal DNA','farther chromosomal DNA'):
  sc.strand(label,[role(label+': sparse covalent tyramide-biotin possible' if stage=='labelled' else label)],label='independent DNA population; no protein/DNA covalent sequence concatenation')
 return sc

def dna_state(label,carrier=False,adapted=False):
 sc=Scene()
 parts=[role(label+' DNA body; original sparse biotin may persist')]
 if adapted:parts=[role('native adapter-derived endA; junction unassigned'),*parts,INDEX_ROLE,role('native adapter-derived endB; junction unassigned')]
 sc.strand('DNA',parts,label=label+' native strand view; opposite end geometry unassigned')
 sc.strand('opposite',[role('opposite genomic/native-copy strand; paired extent unassigned')],label='complementary strand role')
 if carrier:sc.strand('glycogen',[role('glycogen precipitation carrier; downstream partition unassigned')],label='independent carbohydrate polymer; never a copied DNA segment')
 return sc

def clean(g,key,template,label):
 reagent=g.input(key+'_selection',label,MaterialKind.SELECTION)
 kept=g.product(key+'_retained','source recovered '+template.label,MaterialKind.MOLECULAR,state=template.state)
 waste=g.product(key+'_waste','physical removed fraction; class composition unassigned',MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Source physical non-retained fraction; no complete molecular-class depletion claimed.')
 return kept

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input('cells','fixed permeabilized cells; selected SON condition1, other staining conditions are separate variants',MaterialKind.BIOLOGICAL,state=cellular('fixed'),primary=True)
 peroxide=g.input('peroxidase_quench','1.5%H2O2/CMF-PBS1hRT',MaterialKind.CONDITION)
 quenched=g.product('quenched','endogenous-peroxidase-quenched cells',MaterialKind.MOLECULAR,state=cells.state)
 g.react('quench_peroxidase','Quench endogenous peroxidases',StepKind.QUENCHING,((cells,InputRole.TEMPLATE),(peroxide,InputRole.CONDITION)),(quenched,))
 quenched=clean(g,'wash_quench',quenched,'Three PBST centrifugation/pellet washes')
 antibody=g.input('antibody','IgG block, primary20–24h4°C, HRP-secondary20–24h4°C; no-primary control is separate',MaterialKind.SELECTION)
 stained=g.product('stained','antibody-localized HRP population; cells are not selected by staining success',MaterialKind.MOLECULAR,state=cellular('stained'))
 g.react('immunostain','Localize HRP by antibody staining',StepKind.BINDING,((quenched,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(stained,))
 stained=clean(g,'wash_stain',stained,'Three PBST washes after secondary antibody')
 hrp=g.input('hrp','antibody-localized HRP catalyst',MaterialKind.CATALYST)
 tyramide=g.input('tyramide','BufferB:8mM tyramide-biotin stock1:5000+0.003%H2O2; mix equal volume bufferA for final reaction',MaterialKind.MOLECULAR)
 labelled=g.product('labelled','sparsely covalently biotin-labelled chromatin',MaterialKind.MOLECULAR,state=cellular('labelled'))
 g.react('label','Generate localized tyramide radicals10min; no covalent proximity-ligation junction',StepKind.CHEMICAL_CONVERSION,((stained,InputRole.TEMPLATE),(tyramide,InputRole.STRUCTURAL),(hrp,InputRole.CATALYST)),(labelled,))
 labelled=clean(g,'wash_label',labelled,'Three PBST washes after labelling')
 reversal=g.input('reversal','SDS/proteinaseK55°C8–12h,0.2MNaCl65°C20–24h',MaterialKind.CONDITION)
 released=g.product('released','released genomic DNA after deproteinization/reverse-crosslinking',MaterialKind.MOLECULAR,state=dna_state('released'))
 g.react('release','Lyse/proteolyse and reverse crosslinks',StepKind.DEPROTEINIZATION,((labelled,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(released,))
 released=clean(g,'dna_extract',released,'Initialphenol/chloroform extractions; evaporate phenol residues55°C1h')
 rnase=g.input('rnase_a','25µg/mlRNaseA,37°C1h',MaterialKind.CATALYST)
 treated=g.product('rnase_treated','RNase-treated genomic DNA reaction',MaterialKind.MOLECULAR,state=released.state)
 hydrolysate=g.product('rna_hydrolysate','RNA hydrolysis products; physical recovery unassigned',MaterialKind.MOLECULAR)
 g.react('digest_rna','RNaseA treatment after first extraction',StepKind.DIGESTION,((released,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(treated,hydrolysate))
 g.mark_terminal(Terminal(hydrolysate.id,Disposition.UNRESOLVED,'Chemical RNA hydrolysis doesnot itself remove products; downstream extraction/precipitation class recovery unassigned.'))
 released=clean(g,'dna_extract_again',treated,'Re-extract afterRNaseA treatment')
 released=clean(g,'dna_precipitate',released,'SourceNaAc/2volEtOH precipitation−20°C20–120min')
 shear=g.input('sonication','Branson/Bioruptor sonication200–800bp, majority~500bp',MaterialKind.CONDITION)
 fragments=g.product('fragments','random-ended sonicated genomic DNA',MaterialKind.MOLECULAR,state=dna_state('sonicated'))
 g.react('sonicate','Fragment source DNA',StepKind.FRAGMENTATION,((released,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
 inp=g.product('matched_input','10%matched input retained before affinity capture',MaterialKind.MOLECULAR,state=fragments.state)
 capture=g.product('capture_input','remaining90%assigned streptavidin capture',MaterialKind.MOLECULAR,state=fragments.state)
 g.react('split_input','Retain10%matched input; capture remaining DNA',StepKind.SPLIT,((fragments,InputRole.TEMPLATE),),(inp,capture),changes_molecule=False)
 beads=g.input('beads','DynabeadsM270 streptavidin in W&B; selected modified all-at-once captureRT1h/4°C8h',MaterialKind.SELECTION)
 bound=g.product('bound','bead-bound affinity-enriched DNA',MaterialKind.MOLECULAR,state=dna_state('bead-associated'))
 unbound=g.product('unbound','unbound affinity supernatant',MaterialKind.WASTE)
 g.react('capture','Physically capture biotin-bearing DNA',StepKind.AFFINITY_CAPTURE,((capture,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(bound,unbound))
 g.discard(unbound,'Physically removed affinity supernatant; capture is enrichment, not complete label-class separation.')
 bound=clean(g,'wash_capture',bound,'Modified protocol five W&B+0.05%Triton washes55°C; transfer tube each wash')
 release=g.input('capture_release','ProteinaseK+SDS+free biotin overnight55°C',MaterialKind.CONDITION)
 free=g.product('free_pulldown','released pulldown DNA; biotin adduct is not assumed removed',MaterialKind.MOLECULAR,state=dna_state('released pulldown'))
 g.react('release_pulldown','Release captured DNA before extraction',StepKind.DEPROTEINIZATION,((bound,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(free,))
 carrier=g.input('glycogen','0.05µg/µl glycogen precipitation carrier',MaterialKind.MOLECULAR)
 precip=g.input('pulldown_precipitation','Phenol/chloroform; NaAc+2volEtOH12–16h−20°C',MaterialKind.SELECTION)
 recovered=g.product('recovered_pulldown','recovered DNA plus possible independent carrier',MaterialKind.MOLECULAR,state=dna_state('pulldown',True))
 waste=g.product('release_waste','extraction/precipitation non-retained physical fractions',MaterialKind.WASTE)
 g.react('recover_pulldown','Extract and precipitate released pulldown with glycogen',StepKind.CLEANUP,((free,InputRole.TEMPLATE),(carrier,InputRole.SPIKE_IN),(precip,InputRole.SELECTION)),(recovered,waste))
 g.discard(waste,'Source physical fractions; glycogen/DNA class-specific recovery unassigned.')
 for key,template,is_carrier in (('input',inp,False),('pulldown',recovered,True)):
  erat=g.input(key+'_erat','KAPA HyperPrep end-repair/A-tailing reagents; individual enzymes unassigned',MaterialKind.CATALYST)
  repaired=g.product(key+'_repaired','blunt-ended '+key+' DNA',MaterialKind.MOLECULAR,state=dna_state(key+' blunt',is_carrier))
  g.react(key+'_end_repair','Blunt-end source DNA',StepKind.END_REPAIR,((template,InputRole.TEMPLATE),(erat,InputRole.CATALYST)),(repaired,))
  core=Construct([role(key+' end-repaired genomic body; original biotin possible')],name=key+' body')
  da=EndPreparedDuplex(core,five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
  tail_state=dA_tailed_scene(da,label='source3′A;5′status unassigned')
  if is_carrier:tail_state.strand('glycogen',[role('independent glycogen carrier; retention through endprep/ligation unassigned')],label='carbohydrate polymer; never a DNA segment')
  tailed=g.product(key+'_tailed','3′A-tailed '+key+' DNA;5′chemistry unassigned',MaterialKind.MOLECULAR,state=tail_state)
  g.react(key+'_a_tail','Add source3′A tails',StepKind.TAILING,((repaired,InputRole.TEMPLATE),(erat,InputRole.CATALYST)),(tailed,))
  adapters=g.input(key+'_adapters','source6nt-indexed adapters; native strands/member/count/position unassigned',MaterialKind.ADAPTER)
  ligase=g.input(key+'_ligase','KAPA HyperPrep ligation enzymes',MaterialKind.CATALYST)
  adapted=g.product(key+'_adapted','native adapted '+key+' precursor',MaterialKind.MOLECULAR,state=dna_state(key+' adapted',is_carrier,True))
  g.react(key+'_ligate','Ligate source indexed adapters',StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable('Actual3′A substrate known; native chosen adapter strands/junctions not assigned in this source cap.'))
  f=g.input(key+'_forward','native KAPA-compatible PCR primerA',MaterialKind.OLIGO)
  r=g.input(key+'_reverse','native KAPA-compatible PCR primerB',MaterialKind.OLIGO)
  pol=g.input(key+'_polymerase','KAPA HiFi;8–16cycles',MaterialKind.CATALYST)
  lib=g.product(key+'_library',key+' TSA sequencing library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
  original=g.product(key+'_originals','original biotin-bearing DNA/carrier and nonamplified polymers; recovery unassigned',MaterialKind.MOLECULAR,state=adapted.state)
  g.react(key+'_pcr','Amplify adapted '+key+' DNA8–16cycles',StepKind.PCR,((adapted,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(lib,original),geometry=PCRGeometry.unavailable(adapted.state,reason='Actual native adapted precursor; final endpoint is not a template. Original covalent biotin is not automatically copied onto daughter DNA.'))
  g.mark_terminal(Terminal(original.id,Disposition.UNRESOLVED,'Nonamplification does not physically remove originals or carrier; native recovery unassigned.'))
  g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,acquisition='unspecified',limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=('rnase_a','digest_rna','dna_extract_again_selection','dna_extract_again','dna_precipitate_selection','dna_precipitate','cells','peroxidase_quench','quench_peroxidase','wash_quench_selection','wash_quench','antibody','immunostain','wash_stain_selection','wash_stain','hrp','tyramide','label','wash_label_selection','wash_label','reversal','release','dna_extract_selection','dna_extract','sonication','sonicate','split_input','beads','capture','wash_capture_selection','wash_capture','capture_release','release_pulldown','glycogen','pulldown_precipitation','recover_pulldown','input_erat','input_end_repair','input_a_tail','input_adapters','input_ligase','input_ligate','input_forward','input_reverse','input_polymerase','input_pcr','input_library','pulldown_erat','pulldown_end_repair','pulldown_a_tail','pulldown_adapters','pulldown_ligase','pulldown_ligate','pulldown_forward','pulldown_reverse','pulldown_polymerase','pulldown_pcr','pulldown_library')
SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},protocol_graph=protocol_graph,oligos=lambda:())
