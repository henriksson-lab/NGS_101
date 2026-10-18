"""Published CoTECH Table1 component-owned native boundaries; Methods order unverified."""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,InputRole as I,MaterialKind as M,StepKind as K,Terminal,Disposition,SequencingPlan
from pcr import PCRGeometry
import seqprimers as sp
TITLE='CoTECH — source component-qualified joint chromatin/RNA libraries'
NOTES='01_cotech.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-021-01129-z">Xiong et al., Nature Methods2021</a>; published SupplementaryTable1 and SupplementaryNotes. Published defining Methods were not inspected in this capped pass.'
SUMMARY='Published Table1 establishes barcoded pA–Tn5 transfer strands, barcoded oligo(dT), TSO and DNA/RNA library-primer families. Native execution order, shared preamplification and physical branch recovery remain bounded pending the published Methods.'
CAVEAT='INFERRED molecular relationships — this graph is a bounded source-component reconstruction, not a verified execution order. The former preprint-derived pool/redistribute/preamp order is historical and is not active published authority. Verified printed oligos are distinct from unassigned enzyme conditions, loading assignments and first-cycle products. Conditional DNA/RNA branch input classes retain both modalities; their physical allocation is unassigned and primer selectivity does not physically purify modalities. No source absence or validated complete protocol is asserted.'
RTBC=feature('cotech_round1_rt','cell_barcode','whitelist',whitelist='published96RTbarcodes')
T5BC=feature('cotech_round1_t5','cell_barcode','whitelist',whitelist='publishedPA-T5family')
T7BC=feature('cotech_round1_t7','cell_barcode','whitelist',whitelist='publishedPA-T7family')
WELL=feature('cotech_round2','cell_barcode','unknown',note='secondary-well library-primer assignment/width unassigned')
UMI=feature('umi','umi','degenerate',note='printed D8 IUPAC window')
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_INTRO='Published native acquisition was not inspected. First-round barcode/UMI widths are molecular component facts, not an acquired read mask.'
SEQUENCING_ENDING=SEQUENCING_INTRO
READ_SPECS=tuple(sp.read_spec(r,'unassigned','Published run allocation/native read landing not assigned in this bounded pass','undetermined') for r in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
OLIGOS=(('Mosaic-A_primer (MEA)', 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'), ('common annealing primer', '5Phos/CTGTCTCTTATACACATCT'), ('TSO', 'AAGCAGTGGTATCAACGCAGAGTACATrGrG+G'), ('IS primer', 'AAGCAGTGGTATCAACGCAGAGT'), ('P2 primer', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'), ('connector A', 'TCGTCGGCAGCGTCTCCACGC'), ('connector B', 'GTCTCGTGGGCTCGGCTGTCC'), ('Truseq 1st round forward primer', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGTCGGCAGCGTCTCCACGC'), ('Truseq 1st round reverse primer', 'GACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTCTCGTGGGCTCGGCTGTCCC'), ('T5-1 (one of 24 PA-T5 barcode primers)', 'TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG'), ('T7-1 (one of 25 PA-T7 barcode primers)', 'GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG'), ('RT-barcode_01 (one of 96 barcoded oligo dT)', 'TCAGACGTGTGCTCTTCCGATCTAACGTGATDDDDDDDDttttttttttttttttttttttttt'))
T5=OLIGOS[9][1]
T7=OLIGOS[10][1]
RT_ORDERED=OLIGOS[11][1]
TSO='AAGCAGTGGTATCAACGCAGAGTACATGGG'
ME_REV='CTGTCTCTTATACACATCT'
def rt_parts():
 boundary=RT_ORDERED.index('AACGTGAT')
 return [Segment('printed RT5′handle',RT_ORDERED[:boundary]),Segment('RTwell barcode','AACGTGAT',feature=RTBC),Segment('D8 UMI','D'*8,placeholder=True,feature=UMI),Segment('source poly(dT)',RT_ORDERED[boundary+16:].upper())]
def transfer_parts(seq,t5=True):
 barcode='TATAGCCT' if t5 else 'CGAGTAAT';a=seq.index(barcode)
 return [Segment('native adapter prefix',seq[:a]),Segment('source firstwell barcode',barcode,feature=T5BC if t5 else T7BC),Segment('fixed adapter region',seq[a+8:-19]),Segment('mosaic end',seq[-19:])]
def adapter_scene(seq,t5=True):
 sc=Scene();sc.strand('transfer',transfer_parts(seq,t5),label='published transfer strand')
 sc.anneal('nontransfer',[Segment('MErev',ME_REV)],to='transfer',pair=('MErev','mosaic end'),mod5='p',label='source5′phosphorylated common annealing primer');return sc

def mixed_state(stage='initial',indexed=False):
 sc=Scene()
 if stage=='initial':
  sc.strand('DNA',[role('chromatin DNA population')],label='DNA within cell container');sc.strand('RNA',[role('independent poly(A) RNA population; physicalU')],label='RNA within same cell container');return sc
 for key,seq,t5 in (('T5native',T5,True),('T7native',T7,False)):
  sc.strand(key,[*transfer_parts(seq,t5),role('native tagged chromatin strand; gap/other end not reconstructed')]+([role('library-primer index context',molecular_feature=WELL)] if indexed else []),label='INFERRED source-component chromatin reference; '+stage)
 sc.strand('cDNA',[*rt_parts(),role('RNA-complement cDNA; native extent unassigned'),Segment('TSO-complement reference',revcomp(TSO))]+([role('library-primer index context',molecular_feature=WELL)] if indexed else []),label='INFERRED actual ss-cDNA component relationship; native switch/first-cycle scope unverified')
 sc.strand('residual originals',[role('possible original RNA/DNA/nonproductive/background molecules; component recovery unassigned')],label='independent possible original population')
 return sc

def endpoint(name):
 parts=[role('native '+name+' end A'),role('native sourcefirstwell barcode region',length_bp=8,molecular_feature=T5BC if name=='chromatin' else RTBC),role(name+' insert')]
 if name=='chromatin':parts.append(role('opposite firstwell T7 identifier; global order unassigned',length_bp=8,molecular_feature=T7BC))
 else:parts.append(role('UMI identifier; global native readout placement unassigned',length_bp=8,molecular_feature=UMI))
 parts.extend([role('native secondarywell index context',molecular_feature=WELL),role('opposite native end B')])
 return Construct(parts,name='CoTECH '+name+' source-bounded endpoint')
DNA_LIB=endpoint('chromatin');RNA_LIB=endpoint('RNA')
FINAL_LIBRARY=None
FINAL_LIBRARIES=(('Chromatin component reference',DNA_LIB,(),'INFERRED — native global assembly and published reaction order remain unverified.',SEQUENCING_ENDING),('RNA component reference',RNA_LIB,(),'INFERRED — native global assembly and published reaction order remain unverified.',SEQUENCING_ENDING))
def oligos():return tuple(name+' —5′-'+seq+'-3′' for name,seq in OLIGOS)

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input('cells','mixed cell/nucleus DNA and poly(A)RNA; published cell handling unverified',M.BIOLOGICAL,state=mixed_state(),primary=True)
 t5=g.input('t5','published T5-1 transfer/common pair',M.ADAPTER,state=adapter_scene(T5));t7=g.input('t7','published T7-1 transfer/common pair',M.ADAPTER,state=adapter_scene(T7,False))
 rt=g.input('rt_primer','published RTbarcode01, D8 UMI and printed T tract',M.OLIGO,state=Scene());rt.state.strand('RTprimer',rt_parts())
 tso=g.input('tso','published TSO withrGrG+G physical tail',M.OLIGO,state=Scene());tso.state.strand('TSO',[Segment('TSO DNA prefix',TSO[:-3]),Segment('physical TSO terminal GGG',TSO[-3:],note='source rGrG+G; two RNA G and terminal modified +G')])
 native=g.input('firstwell_conditions','native Tn5/RT enzyme choices, temporalorder and conditions unassigned',M.CONDITION)
 first=g.product('firstwell','INFERRED firstwell DNA/ss-cDNA component reference',M.MOLECULAR,state=mixed_state('firstwell sourcecomponent reference'),inferred=True)
 g.react('native_firstwell','INFERRED — native dual labelling boundary; execution order unverified',K.MAGIC,((cells,I.TEMPLATE),(t5,I.STRUCTURAL),(t7,I.STRUCTURAL),(rt,I.STRUCTURAL),(tso,I.STRUCTURAL),(native,I.CONDITION)),(first,),reason='Published Table1 verifies components but defining Methods/reaction order not inspected. No supersededpreprint recipe or enzymes/settings are promoted as published facts.')
 shared=g.input('native_shared_processing','published pooling/redistribution/commonpreamp and recovery specification unread',M.CONDITION)
 prepared=g.product('native_prepared','INFERRED native mixed precursor before branch-specific finishing',M.MOLECULAR,state=mixed_state('native processing boundary'),inferred=True)
 g.react('native_common_boundary','INFERRED — native common processing/recovery boundary',K.MAGIC,((first,I.TEMPLATE),(shared,I.CONDITION)),(prepared,),reason='Published common processing order and actual primer/reagent assignments are unverified; sourcepositive IS/P2/connector components do not establish their execution sequence.')
 aliquot=g.input('aliquot_condition','conditional branch allocation reference; actual physical allocation/amounts/recovery unassigned',M.CONDITION)
 dna=g.product('dna_aliquot','conditional mixed native input class for chromatin finishing',M.MOLECULAR,state=prepared.state);rna=g.product('rna_aliquot','conditional mixed native input class for RNA finishing',M.MOLECULAR,state=prepared.state)
 g.react('branch_inputs','INFERRED — unknown native branch allocation relation; no physical split asserted',K.MAGIC,((prepared,I.TEMPLATE),(aliquot,I.CONDITION)),(dna,rna),reason='Published branch allocation/preparation was not inspected. These are conditional input-class references retaining both modalities, not evidence of a physical aliquot split or modality purification.')
 for key,x,lib in (('dna',dna,DNA_LIB),('rna',rna,RNA_LIB)):
  finish=g.input(key+'_native_conditions','native '+key+' prePCR processing/adaptation/recovery specification unread',M.CONDITION)
  pre=g.product(key+'_precursor','actual source-bounded incoming '+key+' PCR precursor class',M.MOLECULAR,state=mixed_state('unassigned native '+key+' precursor'),inferred=True)
  g.react(key+'_native_finish','INFERRED — native '+key+' substrate preparation boundary',K.MAGIC,((x,I.TEMPLATE),(finish,I.CONDITION)),(pre,),reason='The published precursor-forming operations and their ordered conditions were not inspected. Table1 primer names are not an execution protocol; full actual assembly stays unassigned.')
  f=g.input(key+'_forward','native '+key+' forward-library primer assignment/member unverified',M.OLIGO);r=g.input(key+'_reverse','native '+key+' opposite-library primer assignment/member unverified',M.OLIGO);pol=g.input(key+'_polymerase','native '+key+' library PCRenzyme/cycling unverified',M.CATALYST)
  amp=g.product(key+'_library',lib.name,M.LIBRARY,state=mixed_state('native amplified/cocarried classes',True),inferred=True)
  g.react(key+'_pcr','INFERRED — native '+key+' library amplification relation',K.PCR,((pre,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(pre.state,reason='Incoming native mixed precursor is explicit and source-bounded; completed library is not substituted as PCR template. Exact sourceprimer assignment/native first-cycle geometry not established.'))
  g.sequence(amp,SequencingPlan(library=lib,acquisition='unspecified',limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(('cells','t5','t7','rt_primer','tso','firstwell_conditions','native_firstwell','native_shared_processing','native_common_boundary','aliquot_condition','branch_inputs','dna_native_conditions','dna_native_finish','dna_forward','dna_reverse','dna_polymerase','dna_pcr','dna_library','rna_native_conditions','rna_native_finish','rna_forward','rna_reverse','rna_polymerase','rna_pcr','rna_library'),SOURCE)
 return g.require_publishable()
