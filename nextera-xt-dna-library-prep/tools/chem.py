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

TITLE='Nextera XT DNA library preparation'
NOTES='01_nextera_xt.html'
SOURCE='Illumina Nextera XT Product Documentation15031942v07 April2023; selected standard genomicDNA/XT-v2 indexes/bead normalization. AdapterSequences1000000002694v22 September2025 is later reference only; troubleshooting770-2015-015 qualifies partial adapters and physical ss endpoint.'
SUMMARY='Native partially tagged DNA enters neutralization and the source initial72C extension before indexPCR. Source bead binding, two washes and alkaline strand separation produce a physical single-stranded normalized library.'
CAVEAT+=' Native2023 ATM oligos, gap size, initial-extension geometry and primer strands remain unassigned. Later2025 sequence reference is not proof of native2023 loaded strands. Source XT/v2physical8base indexes are distinguished from acquired cycles and the separate IDT10base option.'
FINAL_LIBRARY=None
SEQUENCING_UNAVAILABLE=''
SOURCE_LIMITATION='Native complete adapter/primer geometry and selected acquisition are unassigned. Source physical XT/v2 indexes have8bases; this does not assign8 acquired cycles. Bead-normalized molecules are single-stranded.'
READ_SPECS=(sp.read_spec('Sequencing output',None,SOURCE_LIMITATION,'include'),)
INDEX_WIDTH=8
PCR_PROGRAM=((72,180),(95,30),('12 cycles',((95,10),(55,30),(72,30))),(72,300))
def partial_state(extended=False):
 sc=Scene()
 for k in ('A','B'):
  sc.strand('native'+k,(role('native partial ATM end '+k+'; loaded strand/ends unassigned'),role('source genomic strand '+k),role('opposite partial boundary; '+('source72C extension applied, completion geometry unassigned' if extended else 'gap/junction geometry unassigned'))),label='actual native tagged DNA strand '+k)
 return sc

def daughter_state(single=False):
 sc=Scene()
 strands=('A',) if single else ('A','B')
 for k in strands:
  sc.strand('daughter'+k,(role('native platform end A; strand arrangement unassigned'),Segment('XT-v2 physical index5 '+k,'N'*INDEX_WIDTH,placeholder=True,feature=I5),role('genomic daughter DNA '+k),Segment('XT-v2 physical index7 '+k,'N'*INDEX_WIDTH,placeholder=True,feature=I7),role('native platform end B; strand arrangement unassigned')),label='normalized single-strand class, polarity unassigned' if single else 'PCR daughter strand '+k)
 sc.note('daughterA','Index feature inventory and physical widths are source-known; drawn order/placement is schematic and unvalidated. Alkaline normalization does not establish which source polarity survives.' if single else 'Two native PCR daughter classes; exact duplex pairing and end ordering are unassigned.')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input('dna','1ng standard genomicDNA; no amplicon-route bead substitution',M.BIOLOGICAL,state=Scene.duplex((Segment('genomicDNA','X'*36,placeholder=True),)),primary=True)
 atm=g.input('atm','ATM tagmentation reagents plus TD;55C5min then10C',M.CATALYST)
 loaded=g.input('loaded_adapters','native ATM-loaded partial adapters; strand composition, ends and geometry unassigned',M.ADAPTER)
 tagged=g.product('tagged','native partial-adapter tagged strands',M.MOLECULAR,state=partial_state())
 nonproductive=g.product('nonproductive','other/nonproductive tagmentation classes; bead recovery unassigned',M.MOLECULAR)
 g.react('tagment','Source solution tagmentation',K.TAGMENTATION,((dna,I.TEMPLATE),(loaded,I.STRUCTURAL),(atm,I.CATALYST)),(tagged,nonproductive))
 g.mark_terminal(Terminal(nonproductive.id,Disposition.UNRESOLVED,'End identity/nonamplifiability does not physically remove a molecule; subsequent cleanup recovery is unassigned.'))
 nt=g.input('nt','5µl NT immediately; room temperature5min',M.CONDITION)
 neutral=g.product('neutralized','NT-treated actual partial native strands',M.MOLECULAR,state=tagged.state)
 g.react('neutralize','Neutralize transposase with NT',K.GENERIC,((tagged,I.TEMPLATE),(nt,I.CONDITION)),(neutral,),changes_molecule=False)
 p5=g.input('p5','selected XT/v2 index primer; native strand/placement unassigned;8base code',M.OLIGO)
 p7=g.input('p7','selected XT/v2 index primer; native strand/placement unassigned;8base code',M.OLIGO)
 npm=g.input('npm','NPM PCR master mix',M.CATALYST)
 extended=g.product('extended','native initial72C-extension material; extent unassigned',M.MOLECULAR,state=partial_state(True))
 g.react('initial_extension','Initial72C3min source PCR-program hold; gap size/completion mechanism unassigned',K.EXTENSION,((neutral,I.TEMPLATE),(p5,I.STRUCTURAL),(p7,I.STRUCTURAL),(npm,I.CATALYST)),(extended,))
 daughter=pcr(g,'index_pcr',extended,daughter_state(),'Source95C30s then12cycles95C10s/55C30s/72C30s and72C5min',p5,p7,npm)
 clean=recover(g,'pcr_cleanup',daughter,'StandardDNA:50µl supernatant+30µlIPB; discard supernatant; two80%EtOH washes; elute52.5µlRSB and transfer50µl')
 aliquot=g.product('normalization_aliquot','20µl cleaned library aliquot',M.MOLECULAR,state=clean.state)
 remainder=g.product('remaining_clean_library','unallocated cleaned library remainder',M.MOLECULAR,state=clean.state)
 g.react('aliquot','Transfer20µl for bead normalization',K.GENERIC,((clean,I.TEMPLATE),),(aliquot,remainder),changes_molecule=False)
 g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Source does not require physical destruction of unused cleaned library.'))
 bound=recover(g,'normalization_bind',aliquot,'Add45µlLNA1/LNB1 mixture,30min; magnet and physically discard unbound supernatant')
 wash1=recover(g,'wash1',bound,'First45µlLNW1 wash; physically discard wash supernatant')
 wash2=recover(g,'wash2',wash1,'Second45µlLNW1 wash; physically discard wash supernatant')
 naoh=g.input('naoh','30µl0.1NNaOH;5min source alkaline elution',M.CONDITION)
 single=g.product('single','alkaline single-stranded normalized library supernatant',M.MOLECULAR,state=daughter_state(True))
 beads=g.product('normalization_beads','remaining normalization beads and retained species; identities unassigned',M.WASTE)
 g.react('alkaline_elution','Alkaline denaturation/elution of source normalized strands',K.STRAND_SEPARATION,((wash2,I.TEMPLATE),(naoh,I.CONDITION)),(single,beads))
 g.discard(beads,'Beads remain behind during source supernatant transfer; species-specific strand partition is unassigned.')
 lns=g.input('lns','fresh plate containing30µlLNS1 per well',M.CONDITION)
 lib=g.product('library','physical single-stranded normalized library',M.LIBRARY,state=single.state)
 g.react('transfer_lns','Transfer30µl alkaline supernatant into LNS1',K.GENERIC,((single,I.TEMPLATE),(lns,I.CONDITION)),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(limitation=SOURCE_LIMITATION))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('dna','atm','loaded_adapters','tagment','nt','neutralize','p5','p7','npm','initial_extension','index_pcr','pcr_cleanup_selection','pcr_cleanup','aliquot','normalization_bind_selection','normalization_bind','wash1_selection','wash1','wash2_selection','wash2','naoh','alkaline_elution','lns','transfer_lns','library')
