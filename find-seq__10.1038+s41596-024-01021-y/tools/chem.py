"""Vendor-revision-owned native molecular declaration."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Disposition, Terminal, SourceDiscrepancy, SourceClaim, AnnotationSlot
from pcr import PCRGeometry
import seqprimers as sp
SEQ_PRIMERS=()
SEQUENCING_INTRO=''
SEQUENCING_ENDING='Illumina-compatible source kit endpoint; complete native primer geometry is unassigned.'
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
def manual_scope(g):
 # This protocol-owned graph explicitly enumerates the inspected native stages.
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
TITLE='FIND-seq — native hydrogel capture and selected-cell transcriptomes'
NOTES='01_find-seq.html'
SOURCE='Shin et al. NatureProtocols2024 doi:10.1038/s41596-024-01021-y; publishedPMC11537836 oligonucleotides and procedure12–102. Selected demonstratedSmart-seq3 route.'
SUMMARY='Acrydite-linked oligo(dT) captures RNA independently of sterically retained genomicDNA. Actual matrix-bound switched cDNA survives detection-PCR sorting; free WTA daughters enter the separately bounded NexteraXT stage.'
CAVEAT='SourceSmart-seq2 alternatives, calibration reactions and bulk/single-cell allocations are distinct from the selected single-cellSmart3 experiment. NativeRNA, genomicDNA, marker-PCR products, matrix originals and possible internal cDNA classes have unresolved recovery unless a physical separation is specified. The referenced NexteraXT revision/kit internals are unread; no canonical index shell is assigned.'
UMI=feature('find_umi','umi','random')
DT=(Segment('matrix spacer','TTTTTT'),Segment('Smart3 DT handle','ACGAGCATCAGCAGCATACGA'),Segment('oligo dT','T'*30),Segment('V anchor','V',placeholder=True))
TSO=(Segment('Smart3 TSO handle','AGAGACAGATTGCGCAATG'),Segment('Smart3 UMI','N'*8,placeholder=True,feature=UMI),Segment('physical terminal riboguanosines','GGG',note='rGrGrG, RNA residues'))
FWD='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG'
REV='ACGAGCATCAGCAGCATACG'
def bead_state():
 sc=population(('genome',(role('sterically retained genomicDNA'),),'independent genome'),('RNA',(role('polyA RNA hybridized to matrix oligo(dT); native association unplaced'),),'independent capturedRNA'),('matrixDT',DT,'actual matrix-linked oligo(dT), not free primer'))
 sc.strands['matrixDT'].mod5='acrydite covalently linked to agarose'
 return sc
def captured(switched=False):
 dna=(*DT,role('first-strand cDNA copied from capturedRNA'))
 if switched:dna=(*dna,*(Segment(p.name,p.top,p.tag,placeholder=p.placeholder,feature=p.feature,note='normalDNA copied from nativeTSO; original RNA/biotin modification not inherited') for p in complement_segments(TSO)))
 sc=population(('genome',(role('sterically retained genomicDNA'),),'independent genome'),('RNA',(role('captured polyA RNA and possible RNA-derived fragments'),),'independent RNA'),('cDNA',dna,'5′matrix-linked actual ss cDNA; not free WTA product'))
 sc.strands['cDNA'].mod5='acrydite covalently linked to agarose'
 return sc
def amplified():
 sc=population(('daughter',(Segment('forward WTA primer-derived end',FWD),Segment('normalDNA UMI copy','N'*8,placeholder=True,feature=UMI),Segment('normalDNA copied terminalTSO riboguanosines','GGG'),role('sense transcript WTA daughter insert'),Segment('copied V anchor','B',placeholder=True),Segment('copied oligo dT','A'*30),role('copied nativeDT primer junction; retained fixedbase context unassigned'),Segment('reverse WTA end',revcomp(REV))),'free double-stranded WTA class; exact end/junction pairing not reconstructed'),('residual',(role('matrix-linked originals, RNA/genome, detection products and other cDNA; partition unresolved'),),'independent residual molecular populations'))
 parts=sc.strands['daughter'].segs
 bottom=complement_segments(parts)
 sc.anneal('daughter complement',bottom,to='daughter',pair=(bottom[-1].name,parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label='free WTA opposite DNA strand; native unknown junctions unpaired')
 return sc
READ_SPECS=(sp.read_spec('Sequencing output',None,'Native selected-cell transcriptomes; source kit/read acquisition remains unassigned','include'),)
def oligos():
 yield oligo('Smart3 matrix oligo(dT)',DT,mods='5′acrydite')
 yield oligo('Smart3 TSO',TSO,mods='5′biotin; terminal rGrGrG')
 yield oligo('Smart3 forward WTA',[Segment('forward',FWD)],mods='terminal A*T*G: two phosphorothioate bonds')
 yield oligo('Smart3 reverse WTA',[Segment('reverse',REV)],mods='source …AC*G*−3′; phosphorothioate markers preserved; terminal typography unresolved')
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cell=g.input('cell','one cell in independent experimental droplet',M.BIOLOGICAL,state=population(('DNA',(role('genomicDNA'),),'genome'),('RNA',(role('polyA RNA'),),'transcriptome')),primary=True)
 dt=g.input('matrix_dt','nativeSmart3 acrydite oligo(dT)',M.OLIGO,state=mol('5′acrydite oligo',*DT))
 agarose=g.input('agarose','allyl agarose/APS/TEMED; source conjugation12–34',M.MOLECULAR)
 matrix=g.product('matrix','covalently primer-functionalized agarose',M.MOLECULAR,state=mol('matrix primer',*DT))
 g.react('conjugate','Conjugate native acrydite oligo to allyl agarose',K.ADAPTER_ATTACHMENT,((dt,I.STRUCTURAL),(agarose,I.TEMPLATE)),(matrix,))
 captured_bead=g.product('captured','coencapsulated lysed genome/RNA bead',M.MOLECULAR,state=bead_state())
 lysis=g.input('lysis','molten agarose encapsulation;LiDS/ProteinaseK55C2h then4C hardening',M.CONDITION)
 g.react('encapsulate_lyse','Capture genome/RNA in functionalized bead',K.GENERIC,((cell,I.TEMPLATE),(matrix,I.STRUCTURAL),(lysis,I.CONDITION)),(captured_bead,),changes_molecule=False)
 x=captured_bead
 for key,label in [('break_emulsion','PFO break; remove oil; cold wash1'),('wash2a','First cold wash2 buffer'),('wash2b','Second cold wash2 buffer'),('rtwash1','First5XRTbuffer wash'),('rtwash2','Second5XRTbuffer wash;100µm filter')]:x=recover(g,key,x,label)
 rt=g.input('rt','MaximaHMinus;roomtemp30min then42C90min',M.CATALYST)
 tso=g.input('tso','nativeSmart3 biotin-TSO withphysicalrGrGrG/N8UMI',M.OLIGO,state=mol('5′biotinTSO',*TSO))
 first=g.product('first','actual matrix-extended cDNA plusRNA/genome',M.MOLECULAR,state=captured())
 g.react('rt_copy','Copy captured RNA from matrix-bound oligo(dT)',K.REVERSE_TRANSCRIPTION,((x,I.TEMPLATE),(matrix,I.STRUCTURAL),(rt,I.CATALYST)),(first,))
 switched=g.product('switched','native switched matrix-bound cDNA mixture',M.MOLECULAR,state=captured(True))
 g.react('switch','Template switch during native RT reaction',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rt,I.CATALYST)),(switched,))
 stop=g.input('stop','EDTA stop;ice1h',M.CONDITION)
 x=op(g,'stop_rt','Stop/cool native RT',K.GENERIC,switched,switched.state,condition=stop)
 for n in range(1,4):x=recover(g,'tween'+str(n),x,'Tween wash'+str(n)+';source centrifugation; final100µm filter')
 markers=g.input('markers','native multiplex detection forward-primer pool; diagnostic species not library handles',M.OLIGO)
 reverse_markers=g.input('marker_reverse','native multiplex detection reverse-primer pool',M.OLIGO)
 probes=g.input('marker_probes','native FAM/VIC hydrolysis probes; diagnostic fluorophores not library identifiers',M.OLIGO)
 taq=g.input('taq','TaqPathProAmp detectionPCR;source55cycles88C denaturation',M.CATALYST)
 measured=g.product('measured','detection-amplified beads plusindependent marker products',M.MOLECULAR,state=captured(True))
 marker_products=g.product('marker_products','independent detectionPCR amplicons/hydrolysed probe products',M.MOLECULAR)
 g.react('marker_pcr','Generate marker fluorescence from actual bead-associated targetDNA/cDNA',K.PCR,((x,I.TEMPLATE),(markers,I.STRUCTURAL),(reverse_markers,I.STRUCTURAL),(probes,I.STRUCTURAL),(taq,I.CATALYST)),(measured,marker_products),geometry=PCRGeometry.unavailable(x.state,reason='Multiplex diagnostic primer/probe geometry unassigned; marker products are not transcriptome library.',symmetric_species=False))
 unresolved(g,marker_products,'Detection products may remain in sorted droplets; their physical recovery through later WTA/beads remains unassigned.')
 gate=g.input('gate','double-positive fluorescence droplet sorting',M.SELECTION)
 selected=g.product('selected','one sorted single-cell droplet',M.MOLECULAR,state=measured.state);reject=g.product('rejected','droplets outside gate',M.WASTE)
 g.react('sort','Physically sort selected single droplets into separate wells',K.SELECTION,((measured,I.TEMPLATE),(gate,I.SELECTION)),(selected,reject));g.discard(reject,'Not deflected into selected collection.')
 cond=g.input('melt','freeze−80C≥2h;thaw;60C10min;aqueous transfer',M.CONDITION)
 x=op(g,'recover_sorted','Melt agarose and transfer aqueous material; matrix linkage retained',K.EXTRACTION,selected,selected.state,condition=cond)
 f=g.input('wta_f','nativeSmart3forward twoPS bonds',M.OLIGO,state=mol('native protectedforward',Segment('forward',FWD,note='terminalA*T*G twoPSbonds')));v=g.input('wta_r','nativeSmart3reverse protectedterminal sourceoligo',M.OLIGO,state=mol('native protectedreverse',Segment('reverse',REV,note='source terminalAC*G* typography;PSmarkers preserved')));pol=g.input('wta_pol','KAPAHiFiReadyMix;empiricallyselected20–26cycle startingrange forscFIND',M.CATALYST)
 x=pcr(g,'wta',x,amplified(),'Amplify actual matrix-linked ss cDNA into free WTA daughters',f,v,pol)
 x=recover(g,'wta_cleanup',x,'0.8XAMPureXP;elute20µl; matrix/original class partition unassigned')
 kit=g.input('xt_kit','referenced NexteraXT kit; exactexecutionrevision unassigned',M.CONDITION)
 library=g.product('library','native NexteraXT-derived selected-cell library',M.LIBRARY,state=mol('native library',role('native platform ends/index selection unassigned'),role('selected-cell daughter cDNA'),role('UMI-bearing tag population where retained',molecular_feature=UMI)))
 g.react('xt_boundary','NexteraXT tagmentation/library preparation per referenced manufacturer instructions',K.MAGIC,((x,I.TEMPLATE),(kit,I.CONDITION)),(library,),reason='ArticleStep102 delegates to manufacturer; exact version-native loaded oligos, neutralization/PCR/cleanup internals are not read in capped source scope. No finalduplex is invented.')
 g.sequence(library,SequencingPlan(limitation=SOURCE_LIMITATION))
 return manual_scope(g)
REVIEWED=('cell','matrix_dt','agarose','conjugate','lysis','encapsulate_lyse','break_emulsion_selection','break_emulsion','wash2a_selection','wash2a','wash2b_selection','wash2b','rtwash1_selection','rtwash1','rtwash2_selection','rtwash2','rt','tso','rt_copy','switch','stop','stop_rt','tween1_selection','tween1','tween2_selection','tween2','tween3_selection','tween3','markers','marker_reverse','marker_probes','taq','marker_pcr','gate','sort','melt','recover_sorted','wta_f','wta_r','wta_pol','wta','wta_cleanup_selection','wta_cleanup','xt_kit','xt_boundary','library')
