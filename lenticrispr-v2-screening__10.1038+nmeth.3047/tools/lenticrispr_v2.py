"""Protocol-owned native declaration; source boundaries stay local."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import crispr as cr
TITLE='lentiCRISPR v2 library construction and validation'
NOTES='01_lenticrispr-v2.html'
SOURCE='Sanjana,Shalem andZhang, <a href="https://doi.org/10.1038/nmeth.3047">NatureMethods2014</a>, published manuscriptPMC4486245 and published MOESM647 SupplementaryMethods286–340; depositedAddgene52961 map for invariantcassette components.'
SUMMARY='The original study builds a pooled plasmid library from synthesized guide oligos, then sequences duplicate barcoded validationPCRs. This is plasmid representation QC; a later screen genomicDNA readout is a separate protocol.'
CAVEAT='Native pooled-library chemistry follows the published supplement. An empty depositedvector map verifies only its invariant hU6/scaffold components; it is not an experimentally cloned guide member. Insert member sequences, exactGibsonjunctions and final QC adapter/barcode strands remain unassigned. Representative clonedguide geometry is INFERRED. A150-cycle kit does not establish the actualRead1/Read2/index allocation. The superseded preprint is not an active source.'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_ENDING='Original source validation: MiSeq withv3 150-cycle kit and independent10%PhiX control; duplicate barcodedPCR populations assess bias. Loaded primers and distribution of150cycles over run roles remain unassigned. The separately sourced compatible screen readout has its own <a href="../joung-one-pcr-crispr-screen-readout__10.1038+nprot.2017.016/joung-one-pcr-crispr-screen-readout.html">Joung one-PCR page</a>.'
ARRAY_F = "TAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"
ARRAY_R = "ACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC"

GUIDE=feature('guide','guide_barcode','unknown')
# Only these invariant hU6/scaffold components are reused; deposited original
# confirms both. The map's filler is not an experimentally cloned guide member.
U6=cr.U6_UPSTREAM+cr.U6_3PRIME
SCAFFOLD=cr.SCAFFOLD_V1

QC_INDEXES={n:feature('qc_replicate_'+str(n),'sample_index','unknown',note='source independently barcoded validation-PCR population; member, width and native strand placement unassigned') for n in (1,2)}

def guide_inventory(sc):
 sc.strand('guide_inventory',[role('source guide-identity inventory; native width/strand placement unassigned',molecular_feature=GUIDE)],label='unplaced identifier inventory; not another physical polymer')

def dna_state(stage):
 sc=Scene()
 if stage=='empty':
  sc.strand('vector',[role('native circular context with cPPT/CTS upstream of hU6'),Segment('hU6 component',U6),role('deposited cloning filler; native cut/end chemistry unassigned'),Segment('scaffold component',SCAFFOLD),role('remaining native vector context')],label='partial empty-vector molecular view; no inserted guide yet')
 elif stage in ('plasmid','genome','screen'):
  prefix=[role('native genomic leftflank; integration sites/copy number unassigned')] if stage!='plasmid' else [role('native circular vector leftcontext; cPPT/CTS upstreamofhU6')]
  suffix=[role('native genomic rightflank')] if stage!='plasmid' else [role('native circular vector rightcontext')]
  parts=prefix+[Segment('hU6 component',U6),Segment('representative extraG','G',inferred=True),Segment('representative20ntguide','N'*20,placeholder=True,inferred=True,feature=GUIDE),Segment('scaffold component',SCAFFOLD)]+suffix
  sc.strand('DNA',parts,label='INFERRED representative inserted-guide locus; integration/circular context not completed')
  sc.strand('other',[role('possible independent unreacted/original DNA, assemblyinsert/filler, oligos, RNA andcellular polymers; class recovery unassigned')])
 else:
  sc.strand('DNA',[role(stage+' native DNA population; full sequence and molecular ends unassigned')]);sc.strand('other',[role('possible independent co-carried DNA/RNA/oligos; class recovery unassigned')])
 if stage.startswith('nativeQC'):guide_inventory(sc)
 return sc

def bounded(g,x,why):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,why))

def oligos():
 yield 'ArrayF:5′-'+ARRAY_F+'-3′; construction primer'
 yield 'ArrayR:5′-'+ARRAY_R+'-3′; construction primer'

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def condition(i,label,kind=M.CONDITION,state=None):return g.input(i,label,kind,state=state)
 def process(i,label,kind,x,state,reagents=(),**kw):
  y=g.product(i+'_out',label,M.MOLECULAR,state=state);g.react(i,label,kind,((x,I.TEMPLATE),*reagents),(y,),**kw);return y
 def cleanup(i,label,x,state=None):
  medium=condition(i+'_selection',label,M.SELECTION);y=g.product(i+'_out','retained '+x.label,M.MOLECULAR,state=state or x.state);w=g.product(i+'_waste','physically unretained fraction',M.WASTE)
  g.react(i,label,K.CLEANUP,((x,I.TEMPLATE),(medium,I.SELECTION)),(y,w),changes_molecule=False);g.discard(w,'Physically unretained sourcefraction; no complete class-specific depletion asserted.');return y
 arraystate=Scene();arraystate.strand('native',[role('microarray-synthesized protected oligonucleotide population; sourcefull-length74nt class',length_bp=74)])
 guide_inventory(arraystate)
 array=g.input('array_oligos','source programmablemicroarray protectedguidepool',M.OLIGO,state=arraystate,primary=True)
 deprotect=condition('deprotection','28–30%ammoniumhydroxide overnight65°C')
 poolstate=Scene();poolstate.strand('pool',[role('released/deprotected single-strand guidepool; sequence/overhangs unassigned',length_bp=74)]);poolstate.strand('other',[role('possible shorter synthesisproducts anddeprotection constituents')])
 guide_inventory(poolstate)
 x=process('release','Cleave anddeprotect sourcearrayproducts',K.CHEMICAL_CONVERSION,array,poolstate,((deprotect,I.CONDITION),))
 dry=condition('dry_te','dry thenresuspend in30µlTE')
 x=process('dry_resuspend','Dry andresuspend synthesisproducts',K.GENERIC,x,x.state,((dry,I.CONDITION),),changes_molecule=False);x=cleanup('oligo_cleanup','SourceQIAquickspincolumn',x)
 f=condition('array_f','ArrayF',M.OLIGO,Construct([Segment('ArrayF',ARRAY_F)]));r=condition('array_r','ArrayR',M.OLIGO,Construct([Segment('ArrayR',ARRAY_R)]));pol=condition('array_pol','PhusionHSFlex; source63°C/15s/20cycles',M.CATALYST)
 ampstate=Scene();ampstate.strand('amplicon',[role('source140bp dsDNApoolPCRclass; native member sequence unassigned',length_bp=140)]);ampstate.strand('other',[role('possible unreacted originals/shorter oligos; subsequent recoveryunassigned')])
 guide_inventory(ampstate)
 amp=g.product('array_pcr_out','source poolPCRdaughters and possible originals',M.MOLECULAR,state=ampstate)
 g.react('array_pcr','Source74ntpool amplification;20cycles',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(x.state,reason='Actual ss synthesispool and possible shortproducts enter PCR; no completed140bp product is used as its template. Native member sequences remain unassigned.'))
 insert=cleanup('insert_gel','Source2%E-GelEX;retain140bpband',amp)
 vector=g.input('vector','native emptylentiCRISPRv2 backbone',M.MOLECULAR,state=dna_state('empty'),primary=True)
 nuc=condition('bsmbi','sourceBsmBI(Fermentas)',M.CATALYST)
 digested=dna_state('BsmBI-opened backbone plus independent excisedfiller')
 digested.strand('filler',[role('independent excisedfiller fragment; stillpresent beforegel recovery')])
 opened=process('digest','Digest sourcevector withBsmBI',K.DIGESTION,vector,digested,((nuc,I.CATALYST),))
 phosph=condition('fastap','sourceFastAP37°C2h',M.CATALYST)
 dephos=Scene();dephos.strand('backbone',[role('native cutvector backbone with both strand5OH; overhangbases unassigned')],mod5='OH');dephos.strand('filler',[role('independent excisedfiller; dephosphorylated mixture, not yet gelremoved')],mod5='OH')
 opened=process('dephosphorylate','Dephosphorylate nativebackbone mixture',K.DEPHOSPHORYLATION,opened,dephos,((phosph,I.CATALYST),))
 retained=Scene();retained.strand('backbone',[role('gelretained nativebackbone, both strand5OH; cutgeometryunassigned')],mod5='OH');retained.strand('possible_carried',[role('possible co-migrating filler/otherDNA; sourcegeldoesnot establishcompletepurity')])
 opened=cleanup('backbone_gel','Source1%E-GelEX;retainbackbone',opened,state=retained)
 gibson=condition('gibson','source20µlNEBGibsonmixture;10nginsert+25ngbackbone',M.CATALYST)
 cloned=g.product('cloned','INFERRED representative pooledguide-bearing vector',M.MOLECULAR,state=dna_state('plasmid'),inferred=True)
 assembly_res=g.product('assembly_other','unassembled/other sourceassemblyclasses',M.MOLECULAR,state=opened.state)
 g.react('gibson_assembly','INFERRED — sourcepooledGibsonassembly',K.MAGIC,((opened,I.TEMPLATE),(insert,I.STRUCTURAL),(gibson,I.CATALYST)),(cloned,assembly_res),reason='Published source specifies pooledGibsonassembly. Individual synthesized74ntmember and nativeoverlap/junction states are not reconstructed; an individualguide-duplex ligation is not substituted.');bounded(g,assembly_res,'Assembly doesnotphysicallyremove unassembledvector/insert/otherproducts; recoveryunassigned.')
 rule=condition('transform_allocation','source0.5µlassembly into25µlelectrocompetentLucigencells')
 aliquot=g.product('transform_aliquot','0.5µlactualassembly aliquot',M.MOLECULAR,state=cloned.state);remainder=g.product('assembly_remainder','unsampledassembly remainder',M.MOLECULAR,state=cloned.state)
 g.react('sample_assembly','Allocate sourceassemblyaliquot',K.SAMPLING,((cloned,I.TEMPLATE),(rule,I.CONDITION)),(aliquot,remainder),changes_molecule=False);bounded(g,remainder,'Sourceunsampledassembly remainder is notprescribeddiscard.')
 selection=condition('bacteria','electrocompetentLucigencells/sourcecarbenicillinplates',M.SELECTION)
 colonies=g.product('colonies','INFERRED source selected colonyplasmidpopulations',M.MOLECULAR,state=dna_state('plasmid'),inferred=True)
 failed=g.product('not_recovered','unrecovered noncolonyfraction',M.WASTE)
 g.react('transform_select','INFERRED — sourceparalleltransformations andcolonyselection',K.SELECTION,((aliquot,I.TEMPLATE),(selection,I.SELECTION)),(colonies,failed),changes_molecule=False);g.discard(failed,'Physicalfraction not recovered withselectedcolonies; no statementthat plasmidfree cells/molecules arechemicallydestroyed.')
 prep=condition('maxiprep_reagent','sourceEndotoxinFreePlasmidMaxiprep(Qiagen)',M.SELECTION)
 library=g.product('plasmid_library','INFERRED source pooledmaxipreppedplasmidlibrary',M.MOLECULAR,state=dna_state('plasmid'),inferred=True);w=g.product('maxiprep_waste','physicallyunretainedextractionfraction',M.WASTE)
 g.react('maxiprep','INFERRED — scrape/combinecolonies andextractsourceplasmids',K.EXTRACTION,((colonies,I.TEMPLATE),(prep,I.SELECTION)),(library,w));g.discard(w,'Sourcephysicalnonretainedextractionfraction; classpurityunassigned.')
 duplicate_rule=condition('duplicate_rule','sourceindependentlybarcoded duplicatevalidationPCRs; aliquotvolumesunassigned')
 aliquots=[g.product('qc_input_'+str(n),'sourcevalidationPCRreplicate'+str(n),M.MOLECULAR,state=library.state) for n in (1,2)]
 rest=g.product('plasmid_rest','plasmidlibrary notallocatedtoQC',M.MOLECULAR,state=library.state)
 g.react('split_qc','Allocate independent sourceQCduplicateinputs',K.SPLIT,((library,I.TEMPLATE),(duplicate_rule,I.CONDITION)),(*aliquots,rest),changes_molecule=False);bounded(g,rest,'Plasmidlibraryremainder outsidevalidationQC; subsequentphysicalfateunassigned.')
 qc_pol=condition('qc_pol','sourcePhusionFlashHighFidelity16cycles',M.CATALYST)
 for n,x in enumerate(aliquots,1):
  f=condition('qc_f_'+str(n),'sourceadapteraddingforwardprimer; strandunassigned',M.OLIGO)
  primer=Scene();primer.strand('native',[role('source native reverse-primer strand; exact sequence unassigned')]);primer.strand('index_inventory',[role('unplaced source technical library-barcode identity',molecular_feature=QC_INDEXES[n])],label='identifier inventory; no width or strand coordinate assigned')
  r=condition('qc_r_'+str(n),'sourcebarcodedreverseprimer; strand/memberunassigned',M.OLIGO,primer)
  daughter=dna_state('nativeQCadapter-bearingdaughterDNA');daughter.strand('index_inventory',[role('unplaced source technical library-barcode identity',molecular_feature=QC_INDEXES[n])],label='identifier inventory; no width or strand coordinate assigned')
  out=g.product('qc_library_'+str(n),'native barcodedQC libraryreplicate'+str(n),M.LIBRARY,state=daughter)
  originals=g.product('qc_originals_'+str(n),'nativeplasmidtemplates andnonselectedclasses',M.MOLECULAR,state=x.state)
  g.react('qc_pcr_'+str(n),'Sourceindependent16cyclevalidationPCR',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(qc_pol,I.CATALYST)),(out,originals),geometry=PCRGeometry.unavailable(x.state,reason='Actual nativeplasmidDNA is theinput. Sourcevalidationprimerstrands/first-cyclegeometry remainunassigned; a compatiblelaterreadout product is not its template.'));bounded(g,originals,'PCR leavesoriginalplasmids/otherclasses chemicallypresent; subsequentphysicalrecoveryunassigned.')
  g.sequence(out,SequencingPlan(limitation='Publishedsource duplicatebarcoded plasmidQC;MiSeqv3 150-cyclekit and10%PhiX. Exactnativeadapter/barcodes, loadedprimers andactualR1/R2/indexcycles remainunassigned. Demultiplexing/trimming is computational, not a physicalmolecularfate.'))
 reviewed=('array_oligos','deprotection','release','dry_te','dry_resuspend','oligo_cleanup_selection','oligo_cleanup','array_f','array_r','array_pol','array_pcr','insert_gel_selection','insert_gel','vector','bsmbi','digest','fastap','dephosphorylate','backbone_gel_selection','backbone_gel','gibson','gibson_assembly','transform_allocation','sample_assembly','bacteria','transform_select','maxiprep_reagent','maxiprep','plasmid_library','duplicate_rule','split_qc','qc_pol','qc_f_1','qc_r_1','qc_pcr_1','qc_f_2','qc_r_2','qc_pcr_2','qc_library_1','qc_library_2')
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()

EXTRA_MOLECULAR_STATES=(('Deposited empty-vector molecular components',Construct(list(dna_state('empty').strands['vector'].segs),name='deposited emptyvector context'),'Source-deposited empty-vector view; hU6/scaffold components and cPPT ordering are source-qualified. Unknown context/filler spans have no assigned length; this is not a cloned-guide member.'),)
