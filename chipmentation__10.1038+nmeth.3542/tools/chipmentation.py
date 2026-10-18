"""Schmidl2015 native bead-bound ChIPmentation; actual tagged DNA before PCR."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,InputRole as I,MaterialKind as M,StepKind as K,Terminal,Disposition,SequencingPlan
import seqprimers as sp
TITLE='ChIPmentation — bead-bound ChIP DNA tagmentation'
NOTES='01_chipmentation.html'
SOURCE='<a href="https://doi.org/10.1038/nmeth.3542">Schmidl et al., Nature Methods2015</a>, Methods ChIPmentation/amplification and published SupplementaryTable2.'
SUMMARY='Antibody-selected bead-bound chromatin is tagmented37°C10min. Physical washes and source ChIP recovery precede native nick-translation/indexPCR; a separate qPCR aliquot chooses enrichment cycles.'
CAVEAT='Selected bead-bound route is distinct from55°C5min purified-ChIP-DNA tagmentation. Native Nextera loaded strands, gaps and exact customATAC primer member/width remain unassigned; no canonical endpoint is substituted. Main Methods explicitly preactivates KAPA98°C45s for successful first-step nick translation. This establishes the operation, not a9nt gap or complete native extension geometry. The24-cycleqPCR assay is not final enrichment cycle count. Source sequencingTable2 contains differing SR/PE runs; no universal index mask inferred.'
SAMPLE=feature('chipmentation_sample','sample_index','unknown',note='Native customNextera index count/member/width unassigned')
FINAL_LIBRARY=Construct([role('native primer-derived end A'),role('ChIP-selected DNA copy'),role('native index-bearing end B',molecular_feature=SAMPLE)],name='source-bounded ChIPmentation daughter library')
SEQ_PRIMERS=()
FINAL_CAPTION='Source-native library end roles; loaded transposon/custom primer sequences not assigned.'
SEQUENCING_INTRO='Published SupplementaryTable2 rows25/26: the same listed K562 ChIPmentation sample has HiSeq2000 SR50 and PE100 runs; other rows remain separate study runs. Index acquisition remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1','50 SR example;100 PE example','target-selected genomic insert; source run variants retained','include'),sp.read_spec('Read 2','100 PE example;not acquired for SR example','opposite genomic insert only in paired-end runs','include'),sp.read_spec('Index 1 (i7)','unassigned','Native acquisition allocation unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native acquisition allocation unassigned','undetermined'))
TAG_MINUTES=10
TAG_TEMP=37
QPCR_CYCLES=24

def chromatin_scene(tagged=False,extended=False):
 sc=Scene()
 for k in ('A','B'):
  p=[role('genomic strand '+k+'; native modifications retained')]
  if tagged:p=[role('native Tn5-transferred end '+k+'; '+('nick translation applied; extent unassigned' if extended else 'gap/nick geometry unassigned')),*p,role('opposite native tagged boundary; composition unassigned')]
  sc.strand(k,p,label='actual native '+('tagged DNA' if tagged else 'chromatin DNA')+' population')
 sc.strand('independent context',[role('protein/antibody/bead and other molecular context; physical recovery unassigned')],label='independent context, not covalent DNA sequence')
 return sc

def daughter_scene():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='normal daughter library DNA',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY))
 sc.strand('residuals',[role('possible original tagged DNA and other co-carried molecules; physical recovery unresolved')],label='independent original/background lineage');return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION):return g.input(i,l,k)
 def op(i,l,k,x,reagents,state):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state);g.react(i,l,k,((x,I.TEMPLATE),*reagents),(y,));return y
 def recover(i,l,x):
  sel=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+' retained fraction',M.MOLECULAR,state=x.state);w=g.product(i+'_waste','physical nonretained fraction',M.WASTE)
  g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));g.discard(w,'Physical source fraction; no complete molecular-class depletion assumed.');return y
 x=g.input('chromatin','crosslinked chromatin',M.BIOLOGICAL,state=chromatin_scene(),primary=True)
 son=inp('sonication','source ChIP sonication');x=op('sonicate','Sonicate crosslinked chromatin',K.FRAGMENTATION,x,((son,I.CONDITION),),x.state)
 ab=inp('antibody','source antibody/capturebeads',M.SELECTION);bound=g.product('bound','bead-bound ChIP selected chromatin',M.MOLECULAR,state=x.state);flow=g.product('flowthrough','IP supernatant',M.WASTE)
 g.react('chip','Immunoprecipitate target-associated chromatin',K.AFFINITY_CAPTURE,((x,I.TEMPLATE),(ab,I.SELECTION)),(bound,flow));g.discard(flow,'Removed ChIP flowthrough; no complete target/background partition assumed.')
 x=recover('pretag_wash','SourceLiClwash followed by two cold10mMTris washes',bound)
 loaded=inp('loaded_ends','native Nextera Tn5-loaded adapter ends; actual strands unassigned',M.ADAPTER);tn5=inp('tn5','Nextera TagmentDNAEnzyme',M.CATALYST)
 x=op('tagment',f'Bead-bound tagmentation{TAG_TEMP}°C{TAG_MINUTES}min',K.TAGMENTATION,x,((loaded,I.STRUCTURAL),(tn5,I.CATALYST)),chromatin_scene(True))
 x=recover('posttag_wash','Two150µl cold sourceWBI/RIPA/TF-WBI washes',x)
 x=recover('final_chip_wash','Chosen ChIP final beadwash; variant-specific buffer unassigned',x)
 elution=inp('chip_elution','source chosen ChIP elution conditions; variant unassigned')
 x=op('elute','Elute tagged ChIP chromatin from beads',K.EXTRACTION,x,((elution,I.CONDITION),),chromatin_scene(True))
 reversal=inp('chip_reversal','source chosen ChIP reversecrosslink/protein handling; exact variant reagents unassigned')
 x=op('reverse_crosslinks','Reverse crosslinks after bead elution',K.DEPROTEINIZATION,x,((reversal,I.CONDITION),),chromatin_scene(True))
 x=recover('dna_purification','Source chosen ChIP DNA purification; native recovery mechanism unassigned',x)
 qamount=inp('q_aliquot','1µl template for10µl qPCR; enrichment remainder separate')
 q=g.product('q_template','qPCR template aliquot',M.MOLECULAR,state=x.state);main=g.product('main_template','remaining native taggedDNA for enrichment',M.MOLECULAR,state=x.state)
 g.react('aliquot','Separate qPCR aliquot and enrichment input',K.PARTITION,((x,I.TEMPLATE),(qamount,I.CONDITION)),(q,main))
 f=inp('forward','customATAC/Nextera forward primer; source member/sequence unassigned',M.OLIGO);r=inp('reverse','customATAC/Nextera reverse primer; source member/sequence unassigned',M.OLIGO)
 pol=inp('kapa','KAPAHiFiHotStartReadyMix preactivated98°C45s',M.CATALYST)
 qex=op('q_nick','qPCR initial72°C5min nick translation',K.EXTENSION,q,((f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),chromatin_scene(True,True))
 qamp=g.product('qamp','24cycle qPCR for cycle estimation',M.QC,state=daughter_scene())
 g.react('qpcr',f'qPCR{QPCR_CYCLES}cycles98/63/72 after98°C30s',K.PCR,((qex,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(qamp,),geometry=PCRGeometry.unavailable(qex.state,reason='Actual native nick-translated aliquot; source custom primer sequence/end layout unassigned.'))
 g.mark_terminal(Terminal(qamp.id,Disposition.QC,'Analytical qPCR cycle choice; not selected sequencing library.'))
 ex=op('main_nick','Enrichment native first-step nick translation; exact fill geometry unassigned',K.EXTENSION,main,((f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),chromatin_scene(True,True))
 choice=inp('cycles','N=rounded-up qPCR Cq; source50µl enrichment reaction')
 amp=g.product('amp','N-cycle enriched library mixture',M.MOLECULAR,state=daughter_scene())
 g.react('enrich','Enrich using qPCR-selected cycle number',K.PCR,((ex,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST),(choice,I.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(ex.state,reason='Actual native extended taggedDNA; no expected finalduplex used as PCR template.'))
 x=recover('size_upper','0.7×AMPure removes long>600bp fraction; retain supernatant',amp)
 x=recover('size_lower','2×AMPure recovery from remaining reaction DNA',x)
 lib=g.product('library',FINAL_LIBRARY.name,M.LIBRARY,state=x.state);g.react('retain','Retain recovered library',K.GENERIC,((x,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(('chromatin','sonication','sonicate','antibody','chip','pretag_wash','pretag_wash_selection','loaded_ends','tn5','tagment','posttag_wash','posttag_wash_selection','final_chip_wash','final_chip_wash_selection','chip_elution','elute','chip_reversal','reverse_crosslinks','dna_purification','dna_purification_selection','q_aliquot','aliquot','forward','reverse','kapa','q_nick','qpcr','qamp','main_nick','cycles','enrich','size_upper','size_upper_selection','size_lower','size_lower_selection','retain','library'),SOURCE)
 return g.require_publishable()
