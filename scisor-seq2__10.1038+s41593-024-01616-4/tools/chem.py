"""Protocol-owned native molecules and source-bounded graph."""
from chemdraw import Construct, Scene, Segment, revcomp, feature, oligo
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE=''
SHOW_READ_CONTENT=False
import seqprimers as sp
def population(*rows):
 sc=Scene()
 for key,parts,label in rows:sc.strand(key,tuple(parts),label=label)
 return sc
def recover(g,key,x,label):
 reagent=g.input(key+'_selection',label,M.SELECTION)
 kept=g.product(key+'_retained','recovered '+x.label,M.MOLECULAR,state=x.state)
 waste=g.product(key+'_waste','physical non-retained supernatant/washes',M.WASTE)
 g.react(key,label,K.CLEANUP,((x,I.TEMPLATE),(reagent,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Actual non-retained fraction; complete species-specific depletion is not asserted.')
 return kept
def unresolved(g,x,reason):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,reason))
def op(g,key,x,label,kind,state,reagent,port=I.CATALYST):
 out=g.product(key+'_out',label,M.MOLECULAR,state=state)
 g.react(key,label,kind,((x,I.TEMPLATE),(reagent,port)),(out,),**({'changes_molecule':False} if kind in (K.GENERIC,K.QUENCHING) else {}))
 return out
def pcr(g,key,x,label,state,f,v,enzyme):
 out=g.product(key+'_out',label,M.MOLECULAR,state=state)
 original=g.product(key+'_original','original/co-carried molecular classes; later physical recovery unassigned',M.MOLECULAR,state=x.state)
 g.react(key,label,K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(v,I.STRUCTURAL),(enzyme,I.CATALYST)),(out,original),geometry=PCRGeometry.unavailable(x.state,reason='Actual incoming native molecular population is supplied; sequential first-cycle geometry is unassigned. Expected final library is not the PCR template.'))
 unresolved(g,original,'Amplification does not erase original templates or physically remove nonamplifying populations; later separation partition remains unassigned.')
 return out
from endprep import EndPreparedDuplex,dA_tailed_scene
TITLE='ScISOr-seq2 — native LAP-CAP cDNA and separate long-read preparations'
NOTES='01_scisor-seq2.html'
SOURCE='Joglekar et al. NatureNeuroscience2024 doi:10.1038/s41593-024-01616-4; publishedPMC11156538 Methods LAP/CAP/PacBio/ONT and originalMOESM1TablesS1–3/MOESM2reportingsummary.'
SUMMARY='Actual native10x cDNA undergoes one-primer asymmetric amplification, then two-primer amplification and physical exome capture. Actual bead-bound DNA/probe mixtures feed on-bead PCR; amplified supernatant is physically recovered before downstream library preparation.'
CAVEAT='Graph independently instantiates PacBio-selected and ONT-selected experiments through LAP-CAP; it does not allocate a single unsupported aliquot to both platforms. Source10xRT/version-native fullstrands, cellbarcode/UMI widths and originalinput endgeometry are unassigned, rather than default16/12. Selectedmouse probe route; human/SnISOr andparallel10xshortread preparations are separate. Exome enrichment/amplification changes representation and does not physically eliminate every unbarcoded/unspliced molecule.'
PARTIAL_READ1='CTACACGACGCTCTTCCGATCT'
PARTIAL_TSO='AAGCAGTGGTATCAACGCAGAGTACAT'
CELL=feature('scisor_cell_barcode','cell_barcode','whitelist',whitelist='source10x3primev3; native members/width unassigned')
UMI=feature('scisor_umi','umi','random',note='native width unassigned in inspected originals')
LIMIT='Native PacBio SMRTbell body/ligatinghairpin ends and precise polymerase circuit geometry remain source-bounded; originalKit2.0 is named but its detailed original is unread. ONTuses nativeSQK-LSK110/ACDE_9110_v110_revC_10Nov2020, not laterLSK114. Longread molecule length is not a fixed acquisitionwidth.'
SEQUENCING_ENDING=LIMIT
READ_SPECS=(sp.read_spec('Sequencing output',None,'PacBioHiFi repeated-template observations or independent ONTPromethION72h experiment; native cell/UMI featurewidths unassigned','include'),)
def cdna(stage='native input',captured=False,probe_present=True):
 sc=population(('sourceA',(role('original native10xDNA strandA; '+stage),role('cellbarcode lineage; width/location unassigned',molecular_feature=CELL),role('UMIlineage; width/location unassigned',molecular_feature=UMI)),'native originalDNA/copy classA'),('sourceB',(role('opposite nativeDNA or carried original templates; '+stage),),'independent original/copy populationB'))
 if 'linear' in stage:sc.strand('asymmetricCopy',(Segment('new PartialRead1-derived5prime',PARTIAL_READ1),role('new preferential ssDNAcopy; native partner/otherend unassigned'),role('copied cellbarcode',molecular_feature=CELL),role('copied UMI',molecular_feature=UMI)),label='new normalDNA one-primer daughter; no symmetric two-primer PCRclaim')
 if 'amplified' in stage:sc.strand('daughter',(Segment('PartialRead1 end',PARTIAL_READ1),role('LAP-CAPdaughter cDNA; heterogeneous splice status'),role('retained cellbarcode',molecular_feature=CELL),role('retained UMI',molecular_feature=UMI),Segment('PartialTSO end reversecomplement',revcomp(PARTIAL_TSO))),label='new normalDNA amplification class; exact duplex geometry bounded')
 if captured:
  if probe_present:sc.strand('probe',(role('nativeexomeprobe species/affinitytag/nativehybrid geometry unassigned'),),label='independent source probe component; bead association begins only at capture')
  sc.strand('blockers',(Segment('PartialRead1 blocker',PARTIAL_READ1),),label='independent native blocking oligo; not covalently joined')
  sc.strand('blocker2',(Segment('PartialTSO blocker',PARTIAL_TSO),),label='independent native blocking oligo')
  sc.note('sourceA','Native cDNA/probe/blocker mixture; bead association occurs at actualcapture operation; native denatured/hybridgeometry unassigned' if probe_present else 'Blocked nativecDNA plusindependent blockers; no exomeprobe orcapturebead added yet')
 return sc

def enriched(g,prefix,source,f,v,pol):
 linear=g.product(prefix+'_linear','native one-primer ssDNA enrichment andoriginals',M.MOLECULAR,state=cdna('linear'))
 g.react(prefix+'_linear_amp','12asymmetric cycles withPartialRead1 only;95C3min;98C20s/64C30s/72C60s',K.EXTENSION,((source,I.TEMPLATE),(f,I.STRUCTURAL),(pol,I.CATALYST)),(linear,))
 x=recover(g,prefix+'_linear_spri',linear,'Source0.8XSPRIselect;two80%EtOHwashes')
 x=pcr(g,prefix+'_exponential',x,'Six two-primer exponentialcycles after actualasymmetric ssDNA mixture',cdna('amplified'),f,v,pol)
 x=recover(g,prefix+'_exponential_spri',x,'Source0.6XSPRIselect;two80%EtOHwashes;30µlEB')
 blockers=g.input(prefix+'_blockers','sourcePartialRead1/PartialTSOblockmix100ng/µl',M.OLIGO,state=(PARTIAL_READ1,PARTIAL_TSO))
 block=g.input(prefix+'_block_conditions','500ngcDNA;95C5min65C5min/hold',M.CONDITION)
 blocked=g.product(prefix+'_blocked','native blocked/denaturedDNA mixture',M.MOLECULAR,state=cdna('amplified',True,False))
 g.react(prefix+'_block','Bind independent blockers toactual enrichedcDNA',K.HYBRIDIZATION,((x,I.TEMPLATE),(blockers,I.STRUCTURAL),(block,I.CONDITION)),(blocked,))
 probes=g.input(prefix+'_probes','sourceMouseAllExon5190-4641/HSQG9611A; native probe assembly unassigned',M.OLIGO)
 hyb=g.input(prefix+'_hyb_conditions','sourceprobe/HybridizationMix/RNaseBlock65C24h',M.CONDITION)
 hybrid=g.product(prefix+'_hybrid','native cDNA/probe/blocker hybrid mixture',M.MOLECULAR,state=cdna('amplified',True))
 g.react(prefix+'_hybridize','Hybridize native exon probes for24h',K.HYBRIDIZATION,((blocked,I.TEMPLATE),(probes,I.STRUCTURAL),(hyb,I.CONDITION)),(hybrid,))
 beads=g.input(prefix+'_beads','M270streptavidin beads prewashed3times;binding30minRT',M.SELECTION)
 captured=g.product(prefix+'_captured','actual on-bead cDNA/probe/blocker population',M.MOLECULAR,state=hybrid.state)
 uncaptured=g.product(prefix+'_uncaptured','physically removed unbound capture supernatant',M.WASTE)
 g.react(prefix+'_capture','Physically capture native probe complexes onM270beads',K.AFFINITY_CAPTURE,((hybrid,I.TEMPLATE),(beads,I.SELECTION)),(captured,uncaptured));g.discard(uncaptured,'Actual unbound supernatant; complete depletion of all unspliced/unbarcoded species is not inferred.')
 x=recover(g,prefix+'_wb1',captured,'Source200µlWB1;RT15min;discardwashsupernatant')
 for n in range(1,4):x=recover(g,prefix+'_wb2_'+str(n),x,'SourceWB2wash'+str(n)+';65C10min each;replacefreshpreheatedWB2')
 water=g.input(prefix+'_resuspend_conditions','sourcewashremoved;resuspendbeads18µlwater',M.CONDITION)
 x=op(g,prefix+'_resuspend',x,'Resuspend actual capturedbeads forPCR',K.GENERIC,x.state,water,I.CONDITION)
 x=pcr(g,prefix+'_capture_pcr',x,'On-bead12cycles98C20s/64C60s/72C3min after95C3min',cdna('amplified'),f,v,pol)
 supernatant=g.product(prefix+'_supernatant','physical amplifiedDNA supernatant fromM270beads',M.MOLECULAR,state=x.state)
 remainder=g.product(prefix+'_bead_remainder','M270beads/probe/template retainedfraction; exact DNApartitionunassigned',M.WASTE)
 select=g.input(prefix+'_supernatant_rule','sourceisolateamplifiedsupernatant fromM270beads',M.SELECTION)
 g.react(prefix+'_release_supernatant','Physically isolate amplifiedsupernatant',K.CLEANUP,((x,I.TEMPLATE),(select,I.SELECTION)),(supernatant,remainder));g.discard(remainder,'Source retains supernatant, leaving actualbeadfraction; no species-specific completeness asserted.')
 return recover(g,prefix+'_capture_spri',supernatant,'Source0.8XSPRIselect aftersupernatant isolation')
def oligos():
 yield oligo('PartialRead1',[Segment('PartialRead1',PARTIAL_READ1)])
 yield oligo('PartialTSO',[Segment('PartialTSO',PARTIAL_TSO)])
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 f=g.input('partial_r1','nativePartialRead1',M.OLIGO,state=PARTIAL_READ1)
 v=g.input('partial_tso','nativePartialTSO',M.OLIGO,state=PARTIAL_TSO)
 pol=g.input('polymerase','sourceKAPAHiFiHotStartReadyMix; exceptseparate10xfront-end',M.CATALYST)
 pac_source=g.input('pac_source','30ngnative10x3primev3cDNA; independentPacBio-directed experimentalpreparation',M.BIOLOGICAL,state=cdna(),primary=True)
 ont_source=g.input('ont_source','30ngnative10x3primev3cDNA; independentONT-directed experimentalpreparation, not same unsupportedaliquot',M.BIOLOGICAL,state=cdna(),primary=True)
 pac=enriched(g,'pb',pac_source,f,v,pol)
 ont=enriched(g,'ont',ont_source,f,v,pol)
 rule=g.input('pac_allocation','source~500ngLAPCAPcDNA forPacBiokit2.0',M.CONDITION)
 pac_aliquot=g.product('pac_aliquot','~500ngnativeenrichedcDNA forPacBio',M.MOLECULAR,state=pac.state)
 pac_remaining=g.product('pac_remaining','PacBioexperiment enrichedDNA notallocated',M.MOLECULAR,state=pac.state)
 g.react('sample_pac','Take actualPacBio-directed sourceamount',K.SAMPLING,((pac,I.TEMPLATE),(rule,I.CONDITION)),(pac_aliquot,pac_remaining),changes_molecule=False);unresolved(g,pac_remaining,'Sampling does not physicallydestroy leftoverCAPcDNA.')
 damage=g.input('damage_mix','SMRTbellExpressTemplatePrepKit2.0;37C30min damage repair',M.CATALYST)
 pac=op(g,'damage_repair',pac_aliquot,'Repair DNA damage onnativeLAPCAPDNA',K.END_REPAIR,pac.state,damage)
 endprep=g.input('end_mix','sourceKit2.0 endrepair/A-tailing20C30min65C30min',M.CATALYST)
 prepared=EndPreparedDuplex(Construct((role('nativebarcodedLAPCAPduplex;source-known3primeA'),)),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
 state=dA_tailed_scene(prepared,label='native3primeA;5primeidentityunassigned')
 state.strand('identifiers',(role('retained cellbarcode',molecular_feature=CELL),role('retained UMI',molecular_feature=UMI)),label='unplacedidentifierfeatures of thisDNA, notfreeoligos')
 pac=op(g,'end_prep',pac,'Sourceendrepair/A-tailing',K.END_REPAIR,state,endprep)
 adapter=g.input('hairpins','nativeSMRTbellExpress2.0 adapters; orderedstrands/ligatingendsunassigned',M.ADAPTER)
 ligase=g.input('pb_ligase','nativeadapterligation20C60min',M.CATALYST)
 closed=cdna('amplified; originalKit2.0SMRTbell closed-template context')
 closed.strand('hairpin_context',(role('two nativehairpin components close oppositeinsertstrands; bodies/orientedjunctionsunassigned'),),label='components ofoneclosedSMRTbell; no invented18nt/A-Tshell')
 library=g.product('pb_ligated','sourceSMRTbell ligationproduct',M.MOLECULAR,state=closed)
 off=g.product('pb_ligation_residual','unligated/partially ligated nativeDNA/offproducts',M.MOLECULAR,state=pac.state)
 g.react('pb_ligate','Ligate nativeKit2.0hairpinadapters',K.LIGATION,((pac,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(library,off),geometry=AdapterLigationGeometry.unavailable('Sourcepaper providesA-tailedinput andSMRTbellKit2.0ligation; exactnativeadapteroverhang/strandgeometryunread, not assumedfrommodern/canonicalkit.'));unresolved(g,off,'Nonclosedeligibility is notphysicaldiscard; actual0.6Xbeadpartition unassigned.')
 library=recover(g,'pb_final_spri',library,'Source0.6XSPRIselect afteradapterligation')
 pblib=g.product('pacbio_library','nativeKit2.0closedSMRTbelllibrary',M.LIBRARY,state=library.state)
 g.react('retain_pacbio','RetainnativePacBiolibrary',K.GENERIC,((library,I.TEMPLATE),),(pblib,),changes_molecule=False)
 g.sequence(pblib,SequencingPlan(limitation=LIMIT+' NativeSMRTbell closure known; exactcircuit/primergeometry unassigned.'))
 rule2=g.input('ont_allocation','source~75fmolLAPCAPcDNA forONT',M.CONDITION)
 ont_aliquot=g.product('ont_aliquot','~75fmol independentONTexperimentLAPCAPDNA',M.MOLECULAR,state=ont.state)
 ont_remaining=g.product('ont_remaining','ONTexperimentCAPDNA notallocated',M.MOLECULAR,state=ont.state)
 g.react('sample_ont','TakeactualsourceONTamount',K.SAMPLING,((ont,I.TEMPLATE),(rule2,I.CONDITION)),(ont_aliquot,ont_remaining),changes_molecule=False);unresolved(g,ont_remaining,'No physicaldestruction impliedby aliquoting.')
 kit=g.input('ont_kit','SQK-LSK110;ACDE_9110_v110_revC_10Nov2020 originalreference',M.CONDITION)
 ontlib=g.product('ont_library','nativeONTmotor-associatedDNA library',M.LIBRARY,state=ont_aliquot.state,inferred=True)
 g.react('ont_boundary','NativeLSK110ampliconsbyligation preparation',K.MAGIC,((ont_aliquot,I.TEMPLATE),(kit,I.CONDITION)),(ontlib,),reason='Article givesexactLSK110guideidentity but orderedmotoradapter/reactionworkups notread in cappedpass. No laterLSK114assembly/geometryorfictionalcovalentproteinpolymer is imported.')
 g.sequence(ontlib,SequencingPlan(acquired_roles=('Nanopore strand',),limitation=LIMIT+' SourcePromethIONFLO-PRO002,72h; Q>7filter is computational, notphysicaldiscard.'))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
# Reviewed operations are explicitly listed for both independent preparations.
_REVIEW_ENRICH=('linear_amp','linear_spri_selection','linear_spri','exponential','exponential_spri_selection','exponential_spri','blockers','block_conditions','block','probes','hyb_conditions','hybridize','beads','capture','wb1_selection','wb1','wb2_1_selection','wb2_1','wb2_2_selection','wb2_2','wb2_3_selection','wb2_3','resuspend_conditions','resuspend','capture_pcr','supernatant_rule','release_supernatant','capture_spri_selection','capture_spri')
REVIEWED=('partial_r1','partial_tso','polymerase','pac_source','ont_source',*(p+'_'+k for p in('pb','ont') for k in _REVIEW_ENRICH),'pac_allocation','sample_pac','damage_mix','damage_repair','end_mix','end_prep','hairpins','pb_ligase','pb_ligate','pb_final_spri_selection','pb_final_spri','retain_pacbio','pacbio_library','ont_allocation','sample_ont','ont_kit','ont_boundary','ont_library')
