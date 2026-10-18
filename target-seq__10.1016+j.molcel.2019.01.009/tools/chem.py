"""Vendor-revision-owned native molecular declaration."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Disposition, Terminal, SourceDiscrepancy, SourceClaim, AnnotationSlot
from pcr import PCRGeometry
import seqprimers as sp
SEQ_PRIMERS=()
SEQUENCING_INTRO=''
SHOW_READ_CONTENT=False
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


from chemdraw import oligo
FINAL_LIBRARY=None
SEQUENCING_UNAVAILABLE=''
SOURCE_LIMITATION='Native selected kit/version/end geometry and complete sequencing-primer landing remain unassigned within this bounded pass.'
def population(*rows):
 sc=Scene()
 for key,parts,label in rows:sc.strand(key,tuple(parts),label=label)
 return sc
def mol(label,*parts):return population(('molecule',parts,label))
def op(g,key,label,kind,x,state,catalyst=None,condition=None):
 inputs=[(x,I.TEMPLATE)]
 if catalyst:inputs.append((catalyst,I.CATALYST))
 if condition:inputs.append((condition,I.CONDITION))
 out=g.product(key+'_out',label,M.MOLECULAR,state=state)
 g.react(key,label,kind,inputs,(out,),**({'changes_molecule':False} if kind==K.GENERIC else {}))
 return out
def unresolved(g,x,reason):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,reason))
def finish(g,x,key,label,reads):
 lib=g.product(key,label,M.LIBRARY,state=x.state)
 g.react(key+'_retain','Retain '+label,K.GENERIC,((x,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(limitation=SOURCE_LIMITATION,**reads))
 return lib
def manual_scope(g):
 # This protocol-owned graph explicitly enumerates the inspected native stages.
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
SOURCE_OLIGO_EXCEPTIONS={'dt':'Native DT is exactly the lowercase IUPAC-VN string in publishedKeyResources; whole-strand source transcription test verifies it. Sourceguard does not match this bundle extraction; no absence inference.'}
TITLE='TARGET-seq — native full-length transcriptome and genotype allocations'
NOTES='01_target-seq.html'
SOURCE='Rodriguez-Meira et al. MolecularCell2019 doi:10.1016/j.molcel.2019.01.009; publishedPMC6436961 STARMethods/keyresources; originalmmc1 full-length96well protocol; mmc3TableS2.'
SUMMARY='Actual genomicDNA, RNA and separately primed first-strand cDNA coexist during native RT/preamp. Physical mixed aliquots feed independent transcriptome and genotype preparations; these are not purified molecule-type partitions.'
CAVEAT='Published Nextera kit statements differ: KeyResources names XTIndexKitv2SetA FC-131-2001/SetC FC-131-2003; detailed full-lengthStep25 names FC-131-1001. No kit/member width is selected from that difference. Selected full-length96well HSPC route. Separate3′TARGET and K562/cellline programs are not merged. STARMethods gives15µl transcriptome aliquot while detailedprotocolStep19 says14µl; the graph selects14µl detailedprocedure and preserves that published discrepancy. Native AccessArray ordered primer ends/barcode side and Nextera index members/widths remain unassigned. RNA/DNA originals, targeted copies and nonproductive products have unresolved physical partition through source recoveries.'
DT_HEAD='AAGCAGTGGTATCAACGCAGAGTAC'
DT=(Segment('OligodT ISPCR handle',DT_HEAD),Segment('oligo dT','T'*30),Segment('VNanchor','VN',placeholder=True))
TSO_BASES='AAGCAGTGGTATCAACGCAGAGTACATGGG'
ISPCR='AAGCAGTGGTATCAACGCAGAGT'
CS1='ACACTGACGACATGGTTCTACA'
CS2='TACGGTAGCAGAGACTTGGTCT'
CELL=feature('target_cell_barcode','cell_barcode','unknown',note='source10base single-direction AccessArray cell index; side/member unassigned')
SAMPLE=feature('transcriptome_i5','sample_index','unknown')
SAMPLE7=feature('transcriptome_i7','sample_index','unknown')
def mixed(stage='native input'):
 sc=population(('gDNA',(role('original genomicDNA and genomic target copies; '+stage),),'independent genomic population'),('RNA',(role('original endogenousRNA and ERCC; processing/retention unassigned'),),'independent RNA populations'))
 if stage!='native input':
  parts=(*DT,role('normal first-strand cDNA; sourceRNA complement'))
  if 'switched' in stage:parts=(*parts,Segment('normalDNA copiedTSO',revcomp(TSO_BASES)))
  sc.strand('cDNA',parts,label='actual ss oligo(dT)-derived cDNA; '+stage)
  sc.strand('target_cDNA',(role('independent mRNA-target-primer ss cDNA and targeted daughters; native ends sample-dependent'),),label='targeted RNA-derived copies')
 return sc
def preamp_state():
 sc=mixed('switched original remains; amplified copies present')
 sc.strand('amplicons',(role('paired genomic/cDNA target amplicons and transcriptome daughters; exact strand inventory source-target-dependent'),),label='independent source preamp daughter classes')
 return sc
def native_xt(stage):
 return population(('xtDNA',(role('native '+stage+' endA; exact strand/ends unassigned',molecular_feature=SAMPLE),role('transcriptome and co-aliquoted targetDNA daughter population'),role('native opposite endB',molecular_feature=SAMPLE7)),'actual native '+stage),('residual',(role('possible originalRNA/gDNA, ss cDNA/primers/nonproductive products; partition unassigned'),),'independent possible carried classes'))
def genotype_state(indexed=False):
 parts=[role('native AccessArrayPE1 end; geometry unassigned')] if indexed else []
 parts.extend((Segment('CS1tag',CS1),role('native genomic-or-cDNA target; distinct populations not joined'),Segment('CS2tag reversecomplement',revcomp(CS2))))
 if indexed:parts.append(role('native PE2/index-bearingprimer inventory;10base barcode side unassigned',molecular_feature=CELL))
 return population(('amplicon',parts,'expected native CS1/CS2-target daughter; physical ends unassigned'),('residual',(role('possible original DNA/RNA/ss cDNA and prior amplicons; recovery unresolved'),),'independent original/copy populations'))
SEQUENCING_ENDING='Separate source libraries: transcriptomeNextSeq75ntsingle-end; genotypeMiSeq151ntR1/10ntindex/151ntR2 with customLNA-sequencingprimers. Physicalbarcode side/complete native sequencing geometry remains unassigned.'
READ_SPECS=(sp.read_spec('Read 1','75 transcriptome;151 genotype','Separate libraries: full-length transcriptome single-end75; customCS1/CS2 target read151','include'),sp.read_spec('Read 2','151 genotype','Target library only; transcriptomeRead2 not acquired','include'),sp.read_spec('Index 1 (i7)','10 genotype; transcriptome unassigned','Native cell barcode read; physical barcode side and primer geometry unassigned','demultiplex'))
def oligos():
 yield oligo('OligodT-ISPCR',DT)
 yield oligo('TSO-LNA',[Segment('TSO',TSO_BASES)],mods='source …ACATrGrG+G;2RNA G then LNA G, not plainDNA reagent')
 yield oligo('ISPCR',[Segment('ISPCR',ISPCR)])
 for label,modified in [('CS1-seq','A+CA+CTG+ACGACATGGTTCTACA'),('CS2-seq','T+AC+GGT+AGCAGAGACTTGGTCT'),('CS1rc-seq','T+GT+AG+AACCATGTCGTCAGTGT'),('CS2rc-seq','A+GAC+CA+AGTCTCTGCTACCGTA')]:
  yield '<b>'+label+'</b>:5′-'+modified+'-3′; source+ denotes LNA; not unmodified PCRtag.'
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cell=g.input('cell','one sortedHSPC; selected full-length96well route',M.BIOLOGICAL,state=mixed(),primary=True)
 ercc=g.input('ercc','independentERCC RNA spike-in',M.MOLECULAR,state=mol('ERCC',role('ERCCRNA')))
 dt=g.input('dt','sourceOligodT-ISPCR in lysisbuffer',M.OLIGO,state=''.join(p.top for p in DT))
 lysis=g.input('lysis','Triton/protease/RNaseinhibitor; sourceRTlysis then72C15min inactivation',M.CONDITION)
 x=g.product('lysate','lysedDNA/RNA/ERCC withnativeDT',M.MOLECULAR,state=mixed())
 g.react('lyse_inactivate','Lyse and heat-inactivate protease before adding RTenzymes',K.GENERIC,((cell,I.TEMPLATE),(ercc,I.TEMPLATE),(dt,I.STRUCTURAL),(lysis,I.CONDITION)),(x,),changes_molecule=False)
 target_rt=g.input('mrna_targets','source sample-specificmRNA-target primerpool TableS2a',M.OLIGO)
 rt=g.input('rt','SMARTScribe42C90min;10×50C2min/42C2min;70C15min',M.CATALYST)
 first=g.product('first','actualss firststrand copies plusnative originals',M.MOLECULAR,state=mixed('firstcopy'))
 g.react('rt_copy','Copy endogenous/ERCCRNA withDT andtargetRTprimers',K.REVERSE_TRANSCRIPTION,((x,I.TEMPLATE),(dt,I.STRUCTURAL),(target_rt,I.STRUCTURAL),(rt,I.CATALYST)),(first,))
 tso=g.input('tso','nativeTSO-LNA …ACATrGrG+G',M.OLIGO,state=mol('physicalTSO',Segment('nativeTSO',TSO_BASES,note='terminal rGrG+G; RNA/RNA/LNA')))
 switched=g.product('switched','actual switchedssDNA and independently targeted material',M.MOLECULAR,state=mixed('switched'))
 g.react('switch','Native template switching during RTprogram',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rt,I.CATALYST)),(switched,))
 isp=g.input('ispcr','sourceISPCR',M.OLIGO,state=ISPCR)
 targets=g.input('targets','sample-specificgDNA/cDNA PCRprimerpairs TableS2a; components explicitly separate fromISPCR',M.OLIGO)
 seqamp=g.input('seqamp','SeqAmp; selectedHSPC22cycles98C15s/67C20s/72C6min',M.CATALYST)
 pre=g.product('preamp_material','native mixed transcriptome/target preamp daughters+originals',M.MOLECULAR,state=preamp_state())
 residual=g.product('preamp_residual','originalRNA/gDNA/ssDNA andnontarget products; recovery unassigned',M.MOLECULAR,state=switched.state)
 g.react('preamp','Amplify actual mixed genomicDNA andss cDNA withISPCR andnative targetpairs',K.PCR,((switched,I.TEMPLATE),(isp,I.STRUCTURAL),(targets,I.STRUCTURAL),(seqamp,I.CATALYST)),(pre,residual),geometry=PCRGeometry.unavailable(switched.state,reason='Actual heterogeneous native originalDNA/ss cDNA input; sample-specific pool and sequential firstcycles unassigned.'))
 unresolved(g,residual,'PCR does not erase original templates; partition through source aliquots and later beads remains unassigned.')
 tx=g.product('tx_aliquot','selected14µl mixed preamp aliquot fortranscriptome',M.MOLECULAR,state=pre.state)
 geno=g.product('geno_aliquot','1.5µl unpurified mixed preamp aliquot forgenotyping',M.MOLECULAR,state=pre.state)
 remain=g.product('remaining','remaining unallocated preamp mixture',M.MOLECULAR,state=pre.state)
 g.react('allocate','Physically take mixed aliquots; no molecular-type purification',K.SPLIT,((pre,I.TEMPLATE),),(tx,geno,remain));unresolved(g,remain,'Remaining mix may be stored; no physical destruction assigned.')
 tx=recover(g,'tx_cleanup',tx,'DetailedStep17–23:14µl+11µlwater+16µlbeads labelled0.6X;twoEtOH;elute8µlEB transfer7.5µl')
 adapters=g.input('xt_adapters','nativeATM-loaded adapters; ordered strands unassigned',M.ADAPTER)
 tn5=g.input('xt_tn5','source700pg DNA/ATM;55C6min;quartervolume',M.CATALYST)
 tagged=g.product('tagged','actualpartial native taggedDNA',M.MOLECULAR,state=native_xt('partially tagged strands'))
 other=g.product('nonproductive','otherTn5 molecular products',M.MOLECULAR)
 g.react('tagment','Install native partialATM adapter ends',K.TAGMENTATION,((tx,I.TEMPLATE),(adapters,I.STRUCTURAL),(tn5,I.CATALYST)),(tagged,other));unresolved(g,other,'Nonproductive ends do not physically remove molecules; bead recovery unassigned.')
 nt=g.input('nt','1.25µlNT neutralization',M.CONDITION)
 tagged=op(g,'neutralize','SourceNT treatment',K.GENERIC,tagged,tagged.state,condition=nt)
 npm=g.input('npm','NPM;sourceinitial72C3min',M.CATALYST)
 extended=op(g,'initial_extension','Apply sourceinitial72C3min; gap size/completion geometry unassigned',K.EXTENSION,tagged,native_xt('initialextension material'),catalyst=npm)
 i5=g.input('i5','nativeNextera i5primer;kitstatement discrepancy/member/widthunassigned',M.OLIGO);i7=g.input('i7','nativeNextera i7primer;member/widthunassigned',M.OLIGO)
 tx=pcr(g,'tx_pcr',extended,native_xt('indexed daughterDNA'),'Source14cycle NPMindexPCR afterinitialextension',i5,i7,npm)
 tx=recover(g,'tx_final_cleanup',tx,'Step30–35 dilute1:1+16µlbeads labelled0.6X;twoEtOH;elute21transfer20µl')
 peers=g.input('tx_peers','independently prepared othercell transcriptome libraries;equalmolarpooling',M.MOLECULAR,state=tx.state)
 pooled=g.product('tx_pooled','physical transcriptome library pool',M.MOLECULAR,state=tx.state)
 g.react('tx_pool','Equimolarpool independently prepared celllibraries',K.POOL,((tx,I.TEMPLATE),(peers,I.TEMPLATE)),(pooled,))
 finish(g,pooled,'transcriptome_library','native full-length transcriptome library',{'acquired_roles':('Read 1',),'read_lengths':{'Read 1':75}})
 p1f=g.input('pcr1f','sourceCS1-tailed target-specificforwardpool TableS2b',M.OLIGO);p1r=g.input('pcr1r','sourceCS2-tailed target-specificreversepool TableS2b',M.OLIGO);kapa=g.input('kapa','KAPA2GRobustHS;20cycles60Cannealing',M.CATALYST)
 nested=pcr(g,'pcr1',geno,genotype_state(),'Diluteunpurifiedaliquot1:2;native20cycleCS1/CS2nestedPCR',p1f,p1r,kapa)
 p2f=g.input('pcr2f','AccessArrayPE1-sidebarcodeprimer;ordered strands unassigned',M.OLIGO);p2r=g.input('pcr2r','AccessArrayPE2-sidebarcodeprimer;barcode side unassigned',M.OLIGO);fast=g.input('fast','FastStartHighFidelity;10cycles60Cannealing',M.CATALYST)
 pcr1_aliquot=g.product('pcr1_aliquot','2.5µlPCR1 physicaltransfer toPCR2',M.MOLECULAR,state=nested.state)
 pcr1_remaining=g.product('pcr1_remaining','unallocatedPCR1remainder',M.MOLECULAR,state=nested.state)
 transfer=g.input('pcr1_transfer','source2.5µlPCR1 aliquot for10µlPCR2',M.CONDITION)
 g.react('sample_pcr1','Take2.5µlPCR1 forPCR2',K.SAMPLING,((nested,I.TEMPLATE),(transfer,I.CONDITION)),(pcr1_aliquot,pcr1_remaining))
 unresolved(g,pcr1_remaining,'SourceunallocatedPCR1 material is not physically destroyed by taking an aliquot.')
 nested=pcr(g,'pcr2',pcr1_aliquot,genotype_state(True),'Transfer2.5µlPCR1 to native10cyclebarcodingPCR2; native10basecellindex',p2f,p2r,fast)
 peers=g.input('genotype_peers','independently barcoded othercell targetlibraries',M.MOLECULAR,state=nested.state)
 pool=g.product('genotype_pool','pooled independently barcoded cellamplicons',M.MOLECULAR,state=nested.state)
 g.react('pool_genotype','Physicalpool onlyafterper-cellbarcodePCR',K.POOL,((nested,I.TEMPLATE),(peers,I.TEMPLATE)),(pool,))
 pool=recover(g,'geno_cleanup',pool,'Sourcepooled0.8XAMPure cleanup')
 finish(g,pool,'genotype_library','native TARGETgenotype library',{'acquired_roles':('Read 1','Read 2','Index 1 (i7)'),'read_lengths':{'Read 1':151,'Read 2':151,'Index 1 (i7)':10}})
 return manual_scope(g)
REVIEWED=('cell','ercc','dt','lysis','lyse_inactivate','mrna_targets','rt','rt_copy','tso','switch','ispcr','targets','seqamp','preamp','allocate','tx_cleanup_selection','tx_cleanup','xt_adapters','xt_tn5','tagment','nt','neutralize','npm','initial_extension','i5','i7','tx_pcr','tx_final_cleanup_selection','tx_final_cleanup','tx_peers','tx_pool','transcriptome_library_retain','transcriptome_library','pcr1f','pcr1r','kapa','pcr1','pcr2f','pcr2r','fast','pcr1_transfer','sample_pcr1','pcr2','genotype_peers','pool_genotype','geno_cleanup_selection','geno_cleanup','genotype_library_retain','genotype_library')
