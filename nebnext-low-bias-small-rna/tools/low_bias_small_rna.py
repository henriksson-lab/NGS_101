"""Vendor-revision-owned native molecular declaration."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Disposition, Terminal, SourceDiscrepancy, SourceClaim, AnnotationSlot
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex,dA_tailed_scene
import seqprimers as sp
SEQ_PRIMERS=()
ACQUIRED_ROLES=()
SEQUENCING_INTRO=''
SEQUENCING_UNAVAILABLE='Native complete primer/adapter strands and selected run are unassigned in this bounded source pass.'
SEQUENCING_ENDING='Illumina-compatible source kit endpoint; complete native primer geometry is unassigned.'
CAVEAT='Unassigned native adapter/primer regions are role tokens. Trim motifs constrain read-through sequence; they do not determine complete kit strands, ends or index placement. Recommendations are distinguished from an acquired run.'
SHOW_READ_CONTENT=False
I5=feature('sample_i5','sample_index','unknown')
I7=feature('sample_i7','sample_index','unknown')
def ss(label,parts):return molecule(label,tuple(parts))
def endpoint(label,payload,left,right):return Construct((left,*payload,right),name=label)
def recover(g,key,template,label,state=None,kind=M.MOLECULAR):
 reagent=g.input(key+'_selection',label,M.SELECTION)
 kept=g.product(key+'_retained','recovered '+template.label,kind,state=template.state if state is None else state)
 waste=g.product(key+'_waste','physically removed supernatant/washes or excluded bead fraction',M.WASTE)
 g.react(key,label,K.CLEANUP,((template,I.TEMPLATE),(reagent,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Physical non-retained fraction; complete species-specific depletion is not inferred.')
 return kept

def pcr(g,key,template,state,label,forward,reverse,polymerase,symmetric=False):
 product=g.product(key+'_product',label+' product',M.MOLECULAR,state=state)
 residual=g.product(key+'_residual','original-template and other co-carried molecular classes; recovery unassigned',M.MOLECULAR,state=template.state)
 inputs=[(template,I.TEMPLATE),(forward,I.STRUCTURAL),(polymerase,I.CATALYST)]
 if reverse is not None:inputs.insert(2,(reverse,I.STRUCTURAL))
 g.react(key,label,K.PCR,inputs,(product,residual),geometry=PCRGeometry.unavailable(template.state,reason='Actual incoming native precursor is shown; first-cycle/end geometry remains unvalidated. Completed library is a separate output.',symmetric_species=symmetric))
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'PCR does not chemically erase original templates or physically purify co-carried polymers; later recovery is unassigned.'))
 return product

def sequence(g,current,endpoint):
 lib=g.product('library','source-qualified sequencing library',M.LIBRARY,state=endpoint)
 g.react('retain_library','Retain recovered library for sequencing',K.GENERIC,((current,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(library=endpoint,limitation=SEQUENCING_UNAVAILABLE))
 return lib

TITLE='NEBNext Low-bias Small RNA'
NOTES='01_nebnext-low-bias-small-rna.html'
SOURCE='<a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale3420.pdf">NEB E3420 manual v2.0_5/26</a>; LV UDI manual v2.0_6/26.'
SUMMARY='Selected post-PCR8A route: native randomized-splint adapter ligations, enzymatic free-adapter processing, source-generated RT primer, two-stage RT bead selection and two-stage short-insert library selection.'
MANUAL_TRIM={'Read 1':'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA','Read 2':'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'}
INPUT=ss('eligible small RNA',(role('native smallRNA less than120nt'),))
INPUT.strands['molecule'].mod5='p';INPUT.strands['molecule'].mod3='OH'
LEFT=role('native LV UDI end A; i5 member/physicalwidth/placement unassigned',molecular_feature=I5)
RIGHT=role('native LV UDI end B; i7 member/physicalwidth/placement unassigned',molecular_feature=I7)
FINAL_LIBRARY=endpoint('native Low-bias SmallRNA8A library',(role('copied smallRNA-derived insert; selected short-insert enrichment'),),LEFT,RIGHT)
FINAL_CAPTION='Native short-insert library. Later June2026 LVUDI reference prints8-base indexes; this does not prove an earlier May2026 executed primer lot. Native selected widths remain unassigned.'
READ_SPECS=(sp.read_spec('Sequencing output','recommendation:56','Vendor says single56-base read is typically sufficient; longer/paired-end possible. No executed read/index recipe assigned.','include'),)

def rna_state(adapters=0):
 parts=[role('native smallRNA insert')]
 if adapters>=1:parts.append(role('native3-prime adapter; strand/end chemistry unassigned'))
 if adapters>=2:parts.insert(0,role('native5-prime adapter; randomized-splint junction unassigned'))
 sc=ss('source RNA substrate',parts)
 if not adapters:sc.strands['molecule'].mod3='OH'
 if adapters<2:sc.strands['molecule'].mod5='p'
 return sc

def cdna_state():
 sc=Scene()
 sc.strand('RNA',rna_state(2).strands['molecule'].segs,label='source RNA template class; native connectivity to modified3-adapter priming moiety unassigned')
 sc.strand('cDNA',[role('generated3-adapter-derived RT primer'),role('antisense copy of smallRNA'),role('native copied5-adapter boundary')],label='ssDNA extension class; free/tethered primer topology unassigned; not completed library')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 rna=g.input('rna','eligible5-primeP/3-primeOH smallRNA; enrichedinput route',M.BIOLOGICAL,state=INPUT,primary=True)
 heat=g.input('heat','§1 RNA denaturation; native reaction preparation',M.CONDITION)
 melted=g.product('melted','denatured eligible RNA',M.MOLECULAR,state=INPUT)
 g.react('denature','RNA denaturation §1',K.THERMAL_TREATMENT,((rna,I.TEMPLATE),(heat,I.CONDITION)),(melted,),changes_molecule=False)
 a3=g.input('adapter3','nativeLB3-prime adapter; strands/activation unassigned',M.ADAPTER)
 l3=g.input('ligase3','NEBNextLB3-prime ligase',M.CATALYST)
 three=g.product('three','3-prime-adapter-ligated RNA',M.MOLECULAR,state=rna_state(1))
 free=g.product('free_adapter','excess unligated3-prime adapter',M.MOLECULAR)
 g.react('ligate3','Ligate native3-prime adapter §2',K.LIGATION,((melted,I.TEMPLATE),(a3,I.STRUCTURAL),(l3,I.CATALYST)),(three,free),geometry=AdapterLigationGeometry.unavailable('SourceRNA3-primeOH known; commercial adapter terminal chemistry and complete strands unassigned.'))
 nuc=g.input('adapter_processing','LBAdaptorRemovalEnzymeMix:30C15min/75C5min',M.CATALYST)
 processed=g.product('processed','ligated RNA after enzyme processing',M.MOLECULAR,state=three.state)
 debris=g.product('processed_adapter','enzymatically processed excess adapter products',M.MOLECULAR)
 g.react('process_adapter','Enzymatically process excess adapter §3',K.DIGESTION,((three,I.TEMPLATE),(free,I.TEMPLATE),(nuc,I.CATALYST)),(processed,debris))
 g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,'Enzymatic processing is not physical removal. Class-specific recovery through later two-stage beads remains unassigned.'))
 a5=g.input('adapter5','nativeLB5-prime adapter',M.ADAPTER)
 splint=g.input('splint','native randomized ligation splint; strand unassigned',M.OLIGO)
 l5=g.input('ligase5','LB5-prime ligase/enhancer/ATP:37C30min',M.CATALYST)
 dual=g.product('dual','both adapters on target RNA',M.MOLECULAR,state=rna_state(2))
 generated=g.product('generated_rt_primer','3-adapter-derived native RT priming component; free/tethered topology unassigned',M.OLIGO,state=ss('native priming component; topology unassigned',(role('3-adapter-derived RT priming moiety; chemistry/connectivity unassigned'),)))
 g.react('ligate5','Ligate5-prime adapter and enzymatically generate RT primer §4/Fig1',K.LIGATION,((processed,I.TEMPLATE),(a5,I.STRUCTURAL),(splint,I.STRUCTURAL),(l5,I.CATALYST)),(dual,generated),geometry=AdapterLigationGeometry.unavailable('Source randomized-splint reaction; native adapter/splint strand geometry unassigned. SourceFig1 states simultaneous3-adapter modification creates RT primer.'))
 rt=g.input('rt','LBRTEnzymeMix/buffer/DTT:42C15min/70C10min',M.CATALYST)
 copied=g.product('copied','native RNA/ss-cDNA reaction',M.MOLECULAR,state=cdna_state())
 g.react('rt_copy','Reverse-transcribe from generated native RT primer §5',K.REVERSE_TRANSCRIPTION,((dual,I.TEMPLATE),(generated,I.STRUCTURAL),(rt,I.CATALYST)),(copied,))
 beads1=g.input('rt_beads1','§6.2 53µl1X beads',M.SELECTION)
 short=g.product('rt_supernatant','retained cDNA-containing supernatant',M.MOLECULAR,state=copied.state)
 large=g.product('large_rt','discarded first-bind large-fragment beads',M.WASTE)
 g.react('rt_first_bind','Bind large products; retain supernatant §6.6',K.SELECTION,((copied,I.TEMPLATE),(beads1,I.SELECTION)),(short,large),changes_molecule=False)
 g.discard(large,'Physically discarded first beads carrying unwanted large fragments.')
 selected=recover(g,'rt_second_bind',short,'§6.7–6.17 add53µl1X beads+40µl isopropanol; retain boundcDNA; discard supernatant/washes')
 f=g.input('forward','LVUDI forward; selectedwell/nativeMay2026width/strand unassigned',M.OLIGO)
 r=g.input('reverse','LVUDI reverse; selectedwell/nativeMay2026width/strand unassigned',M.OLIGO)
 pol=g.input('polymerase','UltraIIQ5:source9–22cycles/62Canneal',M.CATALYST)
 amp=pcr(g,'pcr',selected,FINAL_LIBRARY,'UDI amplification §7',f,r,pol)
 b=g.input('post_first_beads','selected§8A 60µl1.2X beads',M.SELECTION)
 sup=g.product('post_supernatant','retained short-library supernatant',M.MOLECULAR,state=amp.state)
 large2=g.product('post_large','discarded large-fragment first beads',M.WASTE)
 g.react('post_first_bind','§8A.6 retain supernatant anddiscard large-fragment beads',K.SELECTION,((amp,I.TEMPLATE),(b,I.SELECTION)),(sup,large2),changes_molecule=False)
 g.discard(large2,'Physically discarded1.2X first-bind beads; no complete class depletion asserted.')
 final=recover(g,'post_second_bind',sup,'§8A.7–17 add40µl0.8Xbeads; discard supernatant/washes; elute15µl short-library product')
 sequence(g,final,FINAL_LIBRARY)
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('rna','heat','denature','adapter3','ligase3','ligate3','adapter_processing','process_adapter','adapter5','splint','ligase5','ligate5','generated_rt_primer','rt','rt_copy','rt_beads1','rt_first_bind','rt_second_bind_selection','rt_second_bind','forward','reverse','polymerase','pcr','post_first_beads','post_first_bind','post_second_bind_selection','post_second_bind','retain_library','library')
