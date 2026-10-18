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

TITLE='NEBNext Ultra II Directional RNA'
NOTES='01_nebnext-ultra-ii-directional-rna.html'
SOURCE='<a href="https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7760_e7765.pdf">NEB E7760/E7765 original manual v7.0_4/25</a>, selectedSection4 purified mRNA/rRNA-depleted RNA.'
SUMMARY='Source random-primed RNA copies are followed by dUTP second-strand synthesis, native endprep/ligation, USER treatment and source bead recoveries. Native RNA/dU originals remain distinct from daughter DNA.'
CAVEAT+=' The dUTP-marked second strand and USER reagents are source-known; USER cleavage/strand-selectivity and native adapter processing are inferred, not a validated complete-destruction mechanism. No second-strand primer pool is invented from a buffer name.'
SAMPLE=feature('directional_library_index','sample_index','unknown',note='Separately selected oligo kit; single/dual count/width/member unassigned')
FINAL_LIBRARY=Construct((role('native chosen oligo-kit end A'),role('daughter cDNA insert; first-strand lineage reference, direction-selection proof unassigned'),role('native opposite/index end; architecture unassigned',molecular_feature=SAMPLE)),name='native directional RNA endpoint')
FINAL_CAPTION='Native separately selected oligo-kit ends; no dual-eight-base shell or acquired read recipe assigned.'
READ_SPECS=(sp.read_spec('Sequencing output',None,'Source directional-library product; native strand-specific read orientation/primer sites and index acquisition remain unassigned','include'),)
FRAG_MINUTES=(8,15)
RT_PROGRAM=((25,10),(42,15),(70,15))
SECOND_STRAND=(16,60)

def reaction_state(stage,adapted=False):
 sc=Scene()
 first=[role('first-strand DNA copied from source RNA; native sequence unassigned')]
 if adapted:first=[role('native first-strand adapter boundary A; continuity unassigned',inferred=stage=='user'),*first,role('native first-strand adapter boundary B; continuity unassigned',inferred=stage=='user')]
 sc.strand('firstDNA',first,label='normal first-strand DNA; source RNA-complement lineage')
 if stage!='first':
  second=[role('physical second-strand DNA made with dUTP Mix; dU positions unassigned')]
  if stage=='user':second=[role('native dU-bearing second strand and possible USER-cleaved fragments; extent/survival unassigned',inferred=True)]
  if adapted:second=[role('native second-strand adapter context A'),*second,role('native second-strand adapter context B')]
  sc.strand('secondDNA',second,label='independent opposite-strand/processed class; not chemically converted to normal daughter DNA')
 sc.strand('RNA',[role('source RNA and possible RNA-derived products; exact retention/processing unassigned')],label='independent RNA population')
 return sc

def tailed_state():
 prep=EndPreparedDuplex(Construct((role('directional duplex body; native sequences unassigned'),)),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
 sc=dA_tailed_scene(prep,label='source3-primeA; native5-prime chemistry unassigned')
 sc.strands['top'].label='first-strand DNA with3-primeA'
 sc.strands['bottom'].label='opposite dUTP-marked second-strand DNA with3-primeA'
 sc.strand('RNA',[role('independent source RNA-derived population; recovery unassigned')],label='RNA is not a covalent DNA segment')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 rna=g.input('rna','DNA-free purified mRNA or rRNA-depleted RNA; Section4',M.BIOLOGICAL,state=ss('source RNA',(role('source RNA'),)),primary=True)
 buffer=g.input('fragment_buffer','FirstStrandSynthesisReactionBuffer plus kitRandomPrimers;94C sourceinput-dependent fragmentation',M.CONDITION)
 fragments=g.product('fragmented','fragmented/primed RNA',M.MOLECULAR,state=rna.state)
 g.react('fragment','Fragment and prime source RNA §4.1',K.FRAGMENTATION,((rna,I.TEMPLATE),(buffer,I.CONDITION)),(fragments,))
 primer=g.input('random_primers','source kit RandomPrimers; length/sequence unassigned',M.OLIGO)
 ssr=g.input('strand_reagent','source StrandSpecificityReagent; composition unassigned',M.CONDITION)
 rt=g.input('rt','FirstStrandSynthesisEnzymeMix;25C10min/42C15min/70C15min',M.CATALYST)
 first=g.product('first','RNA plus actual first-strand cDNA',M.MOLECULAR,state=reaction_state('first'))
 g.react('rt_copy','Synthesize first-strand cDNA §4.2',K.REVERSE_TRANSCRIPTION,((fragments,I.TEMPLATE),(primer,I.STRUCTURAL),(ssr,I.CONDITION),(rt,I.CATALYST)),(first,))
 dutp=g.input('dutp_buffer','SecondStrandSynthesisReactionBuffer with dUTPMix; no primer formulation assignment',M.CONDITION)
 enzyme=g.input('second_enzyme','SecondStrandSynthesisEnzymeMix;16C60min',M.CATALYST)
 second=g.product('second','normal firstDNA and physical dU-bearing secondDNA',M.MOLECULAR,state=reaction_state('second'))
 g.react('second_copy','Synthesize second-strand DNA using source dUTPMix §4.3',K.EXTENSION,((first,I.TEMPLATE),(dutp,I.CONDITION),(enzyme,I.CATALYST)),(second,))
 clean1=recover(g,'ds_cleanup',second,'§4.4 1.8X144µl beads; twoEtOH washes; recover50µl')
 er=g.input('endprep','UltraII EndPrepEnzymeMix/buffer;20C30min then65C30min',M.CATALYST)
 repaired=g.product('repaired','source end-repaired directional duplex',M.MOLECULAR,state=reaction_state('second'))
 g.react('end_repair','End repair in source EndPrep reaction §4.5',K.END_REPAIR,((clean1,I.TEMPLATE),(er,I.CATALYST)),(repaired,))
 tailed=g.product('tailed','source opposed3-primeA cDNA ends',M.MOLECULAR,state=tailed_state())
 g.react('tail','Add source3-prime dA in same reaction',K.TAILING,((repaired,I.TEMPLATE),(er,I.CATALYST)),(tailed,))
 adapter=g.input('adapter','nonindexed NEBNextAdaptor from separately purchased oligo kit; native strands unassigned',M.ADAPTER)
 ligase=g.input('ligase','UltraII ligationmix/enhancer;20C15min',M.CATALYST)
 adapted=g.product('adapted','native adapter-ligated directional DNA/RNA mixture',M.MOLECULAR,state=reaction_state('second',True))
 g.react('ligate','Ligate native adapter §4.6',K.LIGATION,((tailed,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable('Actual source3-primeA substrate is drawn. Native selected oligo-kit strands, ligating overhang and USER-sensitive structure remain unassigned.'))
 user=g.input('user','source3µlUSER in ligationtube;37C15min',M.CATALYST)
 treated=g.product('user_treated','USER-treated native first/second/RNA molecular classes',M.MOLECULAR,state=reaction_state('user',True),inferred=True)
 cleavage=g.product('user_products','possible USER-generated adapter/dU-strand cleavage products',M.MOLECULAR,state=ss('inferred processing products',(role('native cleavage products; exact identities and recovery unassigned',inferred=True),)))
 g.react('user_treat','INFERRED molecular cleavage scope — USER treatment §4.6.5–6',K.DIGESTION,((adapted,I.TEMPLATE),(user,I.CATALYST)),(treated,cleavage))
 g.mark_terminal(Terminal(cleavage.id,Disposition.UNRESOLVED,'Chemical cleavage/nonamplification is not physical removal. Native RNA/dU-original/adapter-product recovery remains unassigned; complete second-strand destruction is not asserted.'))
 clean2=recover(g,'ligation_cleanup',treated,'§4.7 0.9X87µl beads; washes and native product recovery')
 f=g.input('forward','selected oligo-kit UniversalPCR/i5 forward primer; native sequence/index scheme unassigned',M.OLIGO)
 r=g.input('reverse','selected oligo-kit IndexX/i7 opposite primer; native member/width unassigned',M.OLIGO)
 pol=g.input('polymerase','UltraIIQ5; input-dependent6–13cycles',M.CATALYST)
 amplified=pcr(g,'pcr',clean2,FINAL_LIBRARY,'Copy/amplify actual USER-treated native material §4.8',f,r,pol)
 final=recover(g,'pcr_cleanup',amplified,'§4.9 0.9X beads; twoEtOH washes; recover native library')
 sequence(g,final,FINAL_LIBRARY)
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('rna','fragment_buffer','fragment','random_primers','strand_reagent','rt','rt_copy','dutp_buffer','second_enzyme','second_copy','ds_cleanup_selection','ds_cleanup','endprep','end_repair','tail','adapter','ligase','ligate','user','user_treat','ligation_cleanup_selection','ligation_cleanup','forward','reverse','polymerase','pcr','pcr_cleanup_selection','pcr_cleanup','retain_library','library')
