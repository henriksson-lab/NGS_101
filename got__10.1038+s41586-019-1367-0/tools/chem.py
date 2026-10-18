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
def finish(g,x,key,label,reads):
 lib=g.product(key,label,M.LIBRARY,state=x.state)
 g.react(key+'_retain','Retain '+label,K.GENERIC,((x,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(limitation=SOURCE_LIMITATION,**reads))
 return lib
def manual_scope(g):
 # This protocol-owned graph explicitly enumerates the inspected native stages.
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
import illumina
SOURCE_OLIGO_EXCEPTIONS={'locus_primer':'PublishedSupplementTable5 gives selectedzero-stagger handle andCALR target in separate rows; concatenation is the table-defined primer assembly, components independently source-tested.'}
TITLE='GoT — native linear targeted transcript genotyping'
NOTES='01_got.html'
SOURCE='Nam et al. Nature2019 doi:10.1038/s41586-019-1367-0; publishedPMC6782071 Methods GoT/CircularizationGoT; originalSupplementTable5.'
SUMMARY='Selected3′v2 recovered native first-strand cDNA is amplified before an actual aliquot enters linear genotyping. Native transcriptome remainder is a separate experimental allocation; identifier lineage survives target and index PCR.'
CAVEAT='Selected linear3′v2 route, not a combined3′v2/v3/5′ reagent design. Native10x ordered ends/UMI width are unassigned; source recovered first-strand RNA/DNA populations remain distinct from amplified copies. Circularization and Gibson alternatives are source-known separate experiments summarized explicitly below, with disputed handle assignments retained and full native geometry deferred.'
CELL=feature('cell_barcode','cell_barcode','whitelist',whitelist='native10x3primev2 selected chemistry; members unassigned')
UMI=feature('umi','umi','random')
INDEX=feature('library_index','sample_index','unknown')
P5=illumina.P5
RPI='GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA'
STAGGER='CACCCGAGAATTCCA'
CALR='GCAGCAGAGAAACAAATGAAGGAC'
def native_state(stage):
 return population(('cDNA',(role('native10x endA; '+stage),role('cellbarcode lineage',molecular_feature=CELL),role('nativeUMI width unassigned',molecular_feature=UMI),role('native target-bearing cDNA'),role('native10x endB')),'actual '+stage+' DNA class'),('residual',(role('possible originalRNA, originalssDNA, primers and non-target products; partition unresolved'),),'independent carried native populations'))
def targeted(indexed=False):
 parts=((Segment('P5generic',P5) if indexed else role('native SI-PCR-derived end; exactstrand unassigned')),role('native partialRead1 sequence/identifier geometry unassigned'),role('cellbarcode',molecular_feature=CELL),role('UMI',molecular_feature=UMI),role('variant-containing nativecDNA'),Segment('selected zero-base-stagger handle reverse complement',revcomp(STAGGER)))
 if indexed:parts=(*parts,Segment('sample index','N'*8,placeholder=True,feature=INDEX),Segment('P7 reverse complement',revcomp(illumina.P7)))
 return mol('expected linear PCR daughter; junction/pairedstrand geometry unassigned',*parts)
SEQUENCING_ENDING='Selected linear3′v2 GoT: sourceRead1=26,Read2=98or130,index=8; source expression allocation follows unread native10x recipe.'
CAVEAT+=' CircularGoT is a separate native experiment: targetPCR#1 (optional SF3B1 pre-enrichment),1.3XSPRI,USERII37C1h/65C10min,T4ligation16Covernight/70C15min,lambdaexonuclease37C30min/65C20min,1.3XSPRI;PCR#2 repeats those recoveries/conversion/ligation/digestion;PCR#3 linearizes,SPRI thenP5/RPIindexPCR andPE150MiSeq. Gibson variant substitutes complementary handles/KAPAHiFiReadyMix and50C1hGibson assembly for the U-overhang closure. Lambda digestion makes molecular products; laterSPRI is physical recovery. Native circular junction/topology/ordered source-disputed primer assignment is outside selectedlinear graph. Methods assigns PCR2Fw AGGUCAGTCU/Rv AGACUGACCU; Table5 reverses those labels. MethodsPCR1Rv lacks Table5terminalPS (…CCGATC*T). No molecular choice is asserted for the undrawn circular variants.'
READ_SPECS=(sp.read_spec('Read 1',26,'Native3′v2 cellbarcode+UMI; selected native individual widths unassigned','extract'),sp.read_spec('Read 2','98 or130','targeted variant','include'),sp.read_spec('Index 1 (i7)',8,'sourceRPI library index','demultiplex'))
def oligos():
 yield oligo('P5generic',[Segment('sourceP5generic',P5)])
 yield oligo('RPI-X',[Segment('P7',illumina.P7),Segment('sampleindex','N'*8,placeholder=True,feature=INDEX),Segment('RPI3prime',RPI)])
 yield oligo('selected zero-stagger CALR reverse',[Segment('handle',STAGGER),Segment('CALR',CALR)])
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 first=g.input('first_cdna','selected3′v2 native recovered first-strand10x cDNA; prior version-native kit chemistry bounded',M.BIOLOGICAL,state=native_state('recovered single-stranded firstcopy'),primary=True)
 kit=g.input('native_amp_kit','source10x cDNA amplification reagents/primer mix; one extra3′v2 cycle; native recommendedcyclecount unassigned',M.CONDITION)
 amp=g.product('native_amplified','native amplified10x cDNA mixture',M.MOLECULAR,state=native_state('amplified double-strand copies'))
 g.react('native_amp_boundary','Source native10x cDNA amplification with extra3′v2 cycle',K.MAGIC,((first,I.TEMPLATE),(kit,I.CONDITION)),(amp,),reason='Source delegates native10x version-specific mixture, ordered oligos and base program to manufacturer. Named amplification and sourceextra cycle explicit; no false symmetric unknownprimerpool.')
 amp=recover(g,'native_cleanup',amp,'SourceSPRIselect cleanup after cDNA amplification; ratio/partition unassigned')
 genotype=g.product('genotype_aliquot','3µl native3′v2 cDNA genotyping aliquot',M.MOLECULAR,state=amp.state)
 expression=g.product('expression_remainder','independent remaining native10x cDNA expression allocation',M.MOLECULAR,state=amp.state)
 g.react('allocate','Physically allocate genotyping aliquot and transcriptome remainder',K.SPLIT,((amp,I.TEMPLATE),),(genotype,expression))
 native_lib=g.product('expression_library','native10x expression library',M.LIBRARY,state=native_state('source expression-library daughter'))
 g.react('expression_boundary','Remainder follows standard10x library preparation',K.MAGIC,((expression,I.TEMPLATE),(kit,I.CONDITION)),(native_lib,),reason='Native version-specific fragmentation/adapters/cleanup/primers delegated, unread here. This is a separate aliquot from targeted GoT.')
 g.sequence(native_lib,SequencingPlan(limitation=SOURCE_LIMITATION))
 reamp=g.product('reamplified','3′v2 extra3–4cycle copies plusoriginals',M.MOLECULAR,state=native_state('extra amplified copies'))
 kapa=g.input('kapa','KAPAHiFiHotStartReadyMix+10x primer mix;3–4cycles',M.CONDITION)
 g.react('reamp','Native3′v2 aliquot extra3–4cycle amplification',K.MAGIC,((genotype,I.TEMPLATE),(kapa,I.CONDITION)),(reamp,),reason='Source names commercial primer mix but ordered components and firstcycle geometry unassigned; no invented primer species.')
 x=recover(g,'reamp_cleanup',reamp,'Source cleanup after3′v2 aliquot reamplification; exact ratio unassigned')
 f=g.input('si_pcr','nativeSI-PCR forward; actual ordered strand unassigned',M.OLIGO)
 v=g.input('locus_primer','selected zero-stagger CALR reverse',M.OLIGO,state=STAGGER+CALR)
 pol=g.input('target_pol','native target-PCR reagents; CALR10cycles example; expression-dependent calibration',M.CATALYST)
 x=pcr(g,'target_pcr',x,targeted(),'Enrich actual native CALR-bearing cDNA withSI-PCR andtarget reverse',f,v,pol)
 x=recover(g,'target_cleanup',x,'SourceSPRI recovery removes unincorporated-primer fraction; complete depletion notasserted')
 fp=g.input('p5','sourceP5generic',M.OLIGO,state=P5)
 rp=g.input('rpi','sourceRPI-X; selected8baseindex member unassigned',M.OLIGO,state=mol('RPI-X',Segment('P7',illumina.P7),Segment('index','N'*8,placeholder=True,feature=INDEX),Segment('RPI3prime',RPI)))
 x=pcr(g,'index_pcr',x,targeted(True),'Second sourcePCR completesP7/index; cyclecount unassigned',fp,rp,pol)
 finish(g,x,'genotype_library','native linearGoT genotyping library',{'acquired_roles':('Read 1','Read 2','Index 1 (i7)'),'read_lengths':{'Read 1':26,'Index 1 (i7)':8}})
 return manual_scope(g)
REVIEWED=('first_cdna','native_amp_kit','native_amp_boundary','native_cleanup_selection','native_cleanup','allocate','expression_boundary','expression_library','kapa','reamp','reamp_cleanup_selection','reamp_cleanup','si_pcr','locus_primer','target_pol','target_pcr','target_cleanup_selection','target_cleanup','p5','rpi','index_pcr','genotype_library_retain','genotype_library')
