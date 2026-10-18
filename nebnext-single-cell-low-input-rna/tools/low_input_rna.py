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

TITLE='NEBNext Single Cell / Low Input RNA'
NOTES='01_nebnext-single-cell-low-input-rna.html'
SOURCE='<a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale6420.pdf">NEB E6420 manual v6.0_01/24</a>, selectedSection1 cells/NEBNextUSER-dependent adapter route.'
CAVEAT+=' The source RT/TSO bases are exact; copied TSO boundary and first-cycle strand/amplicon extent are inferred and are not validated reaction proofs.'
SUMMARY='Source RT/TSO oligos create ss cDNA for one-primer amplification; two successive cDNA bead recoveries precede UltraIIFS, native adapter ligation, USER and actual prepared-fragment index PCR.'
RT_PRIMER='AAGCAGTGGTATCAACGCAGAGTAC'+'T'*30+'V'
TSO='GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG'
TSO_DNA=TSO[:-6]
CDNA_PCR_PRIMER='AAGCAGTGGTATCAACGCAGAGT'
RT_HANDLE=RT_PRIMER[:25]
DT=(Segment('source RT25nt handle',RT_HANDLE),Segment('source dT30','T'*30),Segment('source V anchor','V',placeholder=True))
TSO_PARTS=(Segment('source TSO DNA',TSO_DNA),Segment('terminalrGrGrG','GGG',note='three RNA guanosines'))
WTA=Construct((Segment('PCR-derived23nt end',CDNA_PCR_PRIMER),Segment('TSO ACAT','ACAT'),Segment('copiedTSO GGG','GGG',inferred=True),role('RNA-derived full-length cDNA'),Segment('anchor complement','B',placeholder=True),Segment('polyA','A'*30),Segment('RT handle complement',revcomp(RT_HANDLE))),name='source-qualified amplified cDNA')
SAMPLE=feature('library_index','sample_index','unknown',note='Separately chosen NEBNext oligo kit; single/dual layout, native widths/members unassigned')
FINAL_LIBRARY=Construct((role('native separately selected oligo-kit end A'),role('fragmented full-length cDNA insert'),role('native indexed end B; single/dual architecture unassigned',molecular_feature=SAMPLE)),name='native E6420 sequencing endpoint')
FINAL_CAPTION='Unassigned separately purchased adapter/index kit. USER-dependent route selected; no canonical dual-eight-base shell substituted.'
READ_SPECS=(sp.read_spec('Sequencing output','recommendation:2×75','Manual typical paired75-base reads are a recommendation, changeable with insert distribution; no executed run/index mask assigned.','include'),)

def native_fragments(stage):
 sc=Scene()
 sc.strand('fragmentA',[role('native '+stage+' end A'),role('cDNA fragment'),role('native '+stage+' end B')],label='native fragment strand; source-specific boundaries unassigned')
 sc.strand('fragmentB',[role('opposite native '+stage+' cDNA fragment strand')],label='independent complementary DNA strand; actual pairing span unassigned')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 rna=g.input('rna','selected single-cell/low-cell input',M.BIOLOGICAL,state=ss('cell RNA',(role('polyadenylated source RNA'),)),primary=True)
 lysis=g.input('lysis','sourceCellLysisBuffer/RNaseInhibitor and cell RNA preparation §1.1–1.3',M.CONDITION)
 rna.state.strand('gDNA',(role('independent source cellular genomicDNA; downstream fate unassigned'),),label='cellular DNA is not a covalent RNA segment')
 prepared_rna=g.product('prepared_rna','source prepared RNA',M.MOLECULAR,state=rna.state)
 g.react('prepare_rna','Lyse/prepare cell RNA',K.GENERIC,((rna,I.TEMPLATE),(lysis,I.CONDITION)),(prepared_rna,),changes_molecule=False)
 primer=g.input('rt_primer','source NEBNextRT primer with T30V',M.OLIGO,state=Construct(DT))
 rt=g.input('rt','NEBNextSingleCellRT enzyme mix',M.CATALYST)
 first_state=ss('first-strand cDNA',(*DT,role('antisense copy of source RNA'),role('native terminal untemplated boundary')))
 first_state.strand('RNA',(role('independent native RNA template'),),label='source RNA retained; no covalent RNA-cDNA join')
 first_state.strand('gDNA',(role('independent cellular genomicDNA'),),label='DNA carryover unassigned')
 first=g.product('first','actual first-strand RNA/cDNA mixture',M.MOLECULAR,state=first_state)
 g.react('rt_copy','Synthesize first-strand cDNA',K.REVERSE_TRANSCRIPTION,((prepared_rna,I.TEMPLATE),(primer,I.STRUCTURAL),(rt,I.CATALYST)),(first,))
 tso=g.input('tso','source TSO with three terminalrG residues',M.OLIGO,state=Construct(TSO_PARTS))
 copied=tuple(Segment(p.name,p.top,p.tag,p.placeholder,True,feature=p.feature,length_bp=p.length_bp) for p in complement_segments(TSO_PARTS))
 switched_state=ss('actual switched ss-cDNA',(*DT,role('antisense copied RNA'),*copied))
 switched_state.strand('RNA',(role('independent source RNA; thermal fate unassigned'),),label='RNA beside cDNA')
 switched_state.strand('gDNA',(role('independent cellular genomicDNA'),),label='DNA carryover unassigned')
 switched=g.product('switched','source-qualified template-switched ss-cDNA',M.MOLECULAR,state=switched_state,inferred=True)
 g.react('switch','INFERRED copied-end geometry — template-switch using the source TSO',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rt,I.CATALYST)),(switched,))
 wprimer=g.input('wta_primer','one printed23nt cDNA PCR primer',M.OLIGO,state=Construct((Segment('cDNA PCR',CDNA_PCR_PRIMER),)))
 wpol=g.input('wta_polymerase','SingleCellcDNA PCR MasterMix; source11–21cycles',M.CATALYST)
 amp=pcr(g,'wta',switched,WTA,'INFERRED copied-end reference — one-primer full-length cDNA amplification',wprimer,None,wpol,True)
 clean1=recover(g,'cdna_first',amp,'§1.6.2–10 first0.6X60µl beads; wash/elute50µl')
 clean2=recover(g,'cdna_rebind',clean1,'§1.6.11–20 add45µl BeadReconstitutionBuffer to elutedcDNA+samebeads; rebind/wash/final30µl')
 fs=g.input('fs','UltraIIFS enzymes/buffer:37C25min then65C30min',M.CATALYST)
 frag=g.product('fragmented','source FS-fragmented cDNA',M.MOLECULAR,state=native_fragments('fragmented'))
 g.react('fragment','Fragment amplified cDNA §1.8',K.FRAGMENTATION,((clean2,I.TEMPLATE),(fs,I.CATALYST)),(frag,))
 repaired=g.product('repaired','source end-repaired cDNA',M.MOLECULAR,state=native_fragments('repaired'))
 g.react('repair','Repair ends in same FS reaction',K.END_REPAIR,((frag,I.TEMPLATE),(fs,I.CATALYST)),(repaired,))
 prep=EndPreparedDuplex(Construct((role('native FS prepared cDNA insert'),)),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
 tailed=g.product('tailed','source3-primeA-tailed cDNA',M.MOLECULAR,state=dA_tailed_scene(prep,label='source3-primeA;5-prime chemistry unassigned'))
 g.react('tail','Add3-prime dA in same reaction',K.TAILING,((repaired,I.TEMPLATE),(fs,I.CATALYST)),(tailed,))
 adapter=g.input('adapter','selected native NEBNext USER-dependent adapter; separately chosen kit/body unassigned',M.ADAPTER)
 ligase=g.input('ligase','UltraII ligation mix/enhancer:20C15min',M.CATALYST)
 ligated=g.product('ligated','native adapter-ligated cDNA before USER',M.MOLECULAR,state=native_fragments('ligated USER-dependent adapter'))
 g.react('ligate','Ligate selected native NEBNext adapter §1.9',K.LIGATION,((tailed,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable('Source3-primeA substrate; exact separately selected adapter strands, end and USER-sensitive internal structure unassigned.'))
 user=g.input('user','3µl USER:37C15min §1.9.6–7',M.CATALYST)
 treated=g.product('user_treated','USER-processed native cDNA adapters',M.MOLECULAR,state=native_fragments('USER-processed adapter'))
 products=g.product('user_products','USER-cleaved adapter-derived fragments; identity/extent unassigned',M.MOLECULAR)
 g.react('user_treat','Process selected NEBNext adapter with USER',K.DIGESTION,((ligated,I.TEMPLATE),(user,I.CATALYST)),(treated,products))
 g.mark_terminal(Terminal(products.id,Disposition.UNRESOLVED,'USER chemistry is not physical removal; exact released products and recovery unassigned. IndexedUMI adapter route withoutUSER is separate and not recommended by this manual.'))
 clean3=recover(g,'ligation_cleanup',treated,'§1.10 0.8X57µl beads; physically discard supernatant/washes')
 f=g.input('forward','native separately selected oligo-kit forward/universal primer',M.OLIGO)
 r=g.input('reverse','native separately selected indexed opposite primer; single/dual layout unassigned',M.OLIGO)
 pol=g.input('index_polymerase','UltraIIQ5; input-dependent PCR, typical8cycles',M.CATALYST)
 indexed=pcr(g,'index_pcr',clean3,FINAL_LIBRARY,'Amplify actual USER-processed native fragments §1.11',f,r,pol)
 clean4=recover(g,'pcr_cleanup',indexed,'§1.12 0.9X beads; physically discard supernatant/washes')
 sequence(g,clean4,FINAL_LIBRARY)
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('rna','lysis','prepare_rna','rt_primer','rt','rt_copy','tso','switch','wta_primer','wta_polymerase','wta','cdna_first_selection','cdna_first','cdna_rebind_selection','cdna_rebind','fs','fragment','repair','tail','adapter','ligase','ligate','user','user_treat','ligation_cleanup_selection','ligation_cleanup','forward','reverse','index_polymerase','index_pcr','pcr_cleanup_selection','pcr_cleanup','retain_library','library')
