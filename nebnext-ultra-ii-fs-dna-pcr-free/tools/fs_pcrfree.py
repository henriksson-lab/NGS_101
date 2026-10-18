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

TITLE='NEBNext Ultra II FS DNA PCR-free'
NOTES='01_nebnext-ultra-ii-fs-pcr-free.html'
SOURCE='<a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale7430_e7435.pdf">NEB E7430/E7435 original v2.0_6/25</a>, selected high-input350bp route; <a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale7395_e7874_e7876_e7878.pdf">UDI-UMI original v4.0_6/26</a> is a later dated reference.'
SUMMARY='Selected100–500ng/350bp route: coupled FS fragmentation/endprep, native full-length indexed-UMI adapter ligation and two source bead recoveries. No preparation PCR closes the forked ends.'
CAVEAT+=' The2025 source positively states retained single-stranded adapter regions. Later2026 UDI-UMI8/12/20 reference does not identify a2025 lot or native arrangement. The physical library retains forked platform tails; no fully paired shell is manufactured for primer landing.'
FINAL_LIBRARY=None
FINAL_CAPTION='Source PCR-free forked library; native platform-arm sequences, exact stem/junction and identifier positions remain unassigned.'
SOURCE_LIMITATION='Native forked library is the graph endpoint. Sequencing-primer landing and executed read/index acquisition are not assigned; later2026 UDI-UMI settings are explicitly reference scope.'
SEQUENCING_UNAVAILABLE=''
SEQUENCING_ENDING=SOURCE_LIMITATION
READ_SPECS=(sp.read_spec('Sequencing output',None,'Native PCR-free library; executed run unassigned. Later2026 reference:12baseUMI,8baseUDI,Index120withUMI/8without andIndex28; those values do not establish earlier lot identity or native2025 positions.','include'),)
UMI=feature('fs_pcrfree_umi','umi','random',note='Source kit contains UMI; native2025 width/position/multiplicity unassigned; later2026reference12bp')
FS_PROGRAMS={350:(10,37,30,65),450:(8,37,30,65)}
SIZE_SELECT_BEADS_UL={350:65,450:58}

def forked_state():
 core=Segment('genomic insert','X'*36,placeholder=True,note='schematic drawing width; source target350bp is not a fixed molecule length')
 sc=Scene()
 sc.strand('top',[role('native unpaired left adapter arm; i5 context/position unassigned',molecular_feature=I5),core,Segment('retained insert3-primeA','A'),role('native unpaired right adapter arm; i7 context/position unassigned',molecular_feature=I7),role('native UMI-associated context; placement/multiplicity unassigned',molecular_feature=UMI)],label='sealed PCR-free strand; platform tails remain unpaired')
 bottom=[role('opposed native unpaired arm; exact stem unassigned'),*complement_segments((core,)),Segment('opposed insert3-primeA','A'),role('other opposed native unpaired arm; exact stem unassigned')]
 sc.anneal('bottom',bottom,'top',pair=(core.name+"'",core.name),unpaired=tuple(s.name for s in bottom if s.name!=core.name+"'"),label='opposed sealed strand; native adapter junctions/stems unassigned')
 sc.junction('top','native unpaired left adapter arm; i5 context/position unassigned',core.name,'native ligation junction; adapter-side bases unassigned')
 sc.junction('top','retained insert3-primeA','native unpaired right adapter arm; i7 context/position unassigned','native ligation junction; adapter-side bases unassigned')
 sc.note('top','Only genomic core pairing is drawn. Source-known single-strand adapter regions stay forks; native paired adapter stems and exact ligation bases remain unassigned.')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input('dna','100–500ng intact genomicDNA; selected350bp target',M.BIOLOGICAL,state=Scene.duplex((Segment('source DNA','X'*36,placeholder=True),)),primary=True)
 enzyme=g.input('fs','sourceFS enzymes/buffer/DTT; selected37C10min then65C30min',M.CATALYST)
 fragmented=g.product('fragmented','source FS-fragmented DNA',M.MOLECULAR,state=Scene.duplex((Segment('fragment','X'*36,placeholder=True),)))
 g.react('fragment','Fragment DNA in source coupled FS reaction',K.FRAGMENTATION,((dna,I.TEMPLATE),(enzyme,I.CATALYST)),(fragmented,))
 repaired=g.product('repaired','source end-repaired fragments',M.MOLECULAR,state=fragmented.state)
 g.react('repair','End repair in same sourceFS reaction',K.END_REPAIR,((fragmented,I.TEMPLATE),(enzyme,I.CATALYST)),(repaired,))
 prep=EndPreparedDuplex(Construct((role('native end-prepared fragment'),)),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
 tailed=g.product('tailed','source opposed3-primeA substrate',M.MOLECULAR,state=dA_tailed_scene(prep,label='source3-primeA; native5-prime chemistry unassigned'))
 g.react('tail','Add3-prime dA in same FS reaction',K.TAILING,((repaired,I.TEMPLATE),(enzyme,I.CATALYST)),(tailed,))
 adapter=g.input('adapter','full-length UDI-UMI adapter;2025 native strands/end/index arrangement unassigned',M.ADAPTER)
 ligase=g.input('ligase','UltraII LigationMasterMix/enhancer;20C15min',M.CATALYST)
 ligated=g.product('ligated','native sealed forked PCR-free DNA',M.MOLECULAR,state=forked_state())
 other=g.product('other_ligation','other/nonlibrary ligation products',M.MOLECULAR)
 g.react('ligate','Ligate native full-length indexedUMI adapters',K.LIGATION,((tailed,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(ligated,other),geometry=AdapterLigationGeometry.unavailable('Known3-primeA substrate is drawn; exact native2025 adapter strands/overhang and pairedstem are unassigned. Later2026 compatible-kit TA description is reference scope, not oldlot proof.'))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Nonproductive molecules are not removed merely by failing cluster/PCR eligibility. Actual recovery through the two source beads remains unassigned.'))
 first=recover(g,'cleanup',ligated,'§1.3.1–13 add40µl0.1XTE then50µlbeads; discard supernatant/washes; elute102µl and transfer100µl')
 final=recover(g,'size_select',first,'selected350bp §1.3.14–24 add65µlbeads to100µl; discard supernatant/washes; elute22µl and transfer20µl')
 library=g.product('library','recovered forked PCR-free library',M.LIBRARY,state=final.state)
 g.react('retain_library','Retain source forked library; qPCR is analytical quality control outside preparation',K.GENERIC,((final,I.TEMPLATE),),(library,),changes_molecule=False)
 g.sequence(library,SequencingPlan(limitation=SOURCE_LIMITATION))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('dna','fs','fragment','repair','tail','adapter','ligase','ligate','cleanup_selection','cleanup','size_select_selection','size_select','retain_library','library')
