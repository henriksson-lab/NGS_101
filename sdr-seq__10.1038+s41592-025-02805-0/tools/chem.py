"""SDR-seq native V2 mixed droplets and source-specific library PCRs."""
from chemdraw import Scene,Segment,feature,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import seqprimers as sp
TITLE='SDR-seq — source V2 matched DNA/RNA libraries'
NOTES='01_sdr_seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-025-02805-0">Lindenhofer etal.,NatureMethods2025</a>; publishedMOESM3Tables1–3 and cited<a href="https://doi.org/10.17504/protocols.io.6qpvr9q43vmk/v1">detailedprotocolv1,August11,2025</a>.'
SUMMARY='In situ barcoded RT precedes Tapestri v2 celllysis and mixed droplet multiplexPCR. The recoveredmixedDNA/cDNA population supplies both R2N-DNA and R2-RNA indexing reactions; aliquoting does not physicallypurify modalities.'
CAVEAT='SelectedV2 source design, not V1proof-of-concept which usedR2N onbothbranches. PublishedTable3 usesV2RTprimer names forPOP althoughTable1 suppliesV1 sequences withoutGACACGTC; actualPOP sequence remains unresolved. Detailedv1 RNAindexPCR table labels reverseprimerR2N whileits narrative and publishedV2Table2 specifyRNA R2; thismodel usesR2. MainMethods RTmix says500mMdNTPs whiledetailedprotocol table states500µM perbase; selecteddetailed500µM. MainMethods namesMaximaHMinus EP0752 whiledetailedprotocolMaterials listsEP0751; enzymefamilyisretainedwithoutresolvingcatalogidentity. NativeTapestribead sequences/linker topology/PCRpolymerase and unreadv2guide workups stayunassigned. Readlengths beloware explicitdetailedv1recommendations, not inferred from physicalindices. No geneediting experiment is partofthislibrary graph.'
CS='GTACTCGCAGTAGTC'
RNA_FORWARD=CS+'GACACGTC'
R2N='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
R2='GTGACTGGAGTTCAGACGTGT'
P5='AATGATACGGCGACCACCGAGATCTACAC'
P7='CAAGCAGAAGACGGCATACGAGAT'
R1N_INDEX_ARM='TCGTCGGCAGCGTC'
R2N_INDEX_ARM='GTCTCGTGGGCTCGG'
RT_SAMPLE=feature('sdr.rt_sample','sample_index','whitelist',whitelist='publishedV2RTsamplebarcodefamily',note='representativeRT_primer1_8 TCGCCTTA')
UMI=feature('sdr.umi','umi','random')
CELL1=feature('sdr.cell1','cell_barcode','unknown',group='sdr.cell',part='1')
CELL2=feature('sdr.cell2','cell_barcode','unknown',group='sdr.cell',part='2')
I5=feature('sdr.i5','sample_index','whitelist',whitelist='sourceV2Nextera_R1N_S50x')
I7=feature('sdr.i7','sample_index','whitelist',whitelist='sourceV2N70x')
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_ENDING='Detailedv1 recommendsDNA161/161/i5=8/i7=8 andRNA76/46/i5=8/i7=8; exactexperimentrunmask notuniversallyassigned. Gene/panelgenotype/qualityfilters arecomputational.'
READ_SPECS=(sp.read_spec('Read 1','161DNA /76RNA','dna:cellbarcode+amplicon;RNA:cellbarcode,RTsamplebarcode,N8UMI','extract_then_include'),sp.read_spec('Read 2','161DNA /46RNA','dna/RNAamplicon','include'),sp.read_spec('Index 1 (i7)',8,'librarysampleindex bothlibraries','demultiplex'),sp.read_spec('Index 2 (i5)',8,'librarysampleindex bothlibraries','demultiplex'))
def rt_parts():return [Segment('V2CS/amplificationhandle',RNA_FORWARD),Segment('sourceRTsamplebarcodeexample','TCGCCTTA',feature=RT_SAMPLE),Segment('N8UMI','N'*8,feature=UMI),Segment('dT16','T'*16),Segment('VNanchor','VN',placeholder=True)]
def cell_parts():return [role('cellBC1; widthunassigned',molecular_feature=CELL1),role('variablelinker14–17bases; selectedmemberunassigned'),role('cellBC2; widthunassigned',molecular_feature=CELL2),role('source15bpconstantlinker; identitywithCSnotestablished',length_bp=15)]
def ss(name,parts):
 sc=Scene();sc.strand(name,parts,label='source5′→3′DNAoligo');return sc

def mixture(stage):
 sc=Scene()
 sc.strand('originalgenome',[role('original genomicDNA; sourcefixation/nativefragmentextentunassigned')],label='independentoriginalDNA')
 sc.strand('originalRNA',[role('originalpolyARNAwithphysicalU orretainedprocessingproducts; RNAfateunassigned')],label='independentRNA')
 if stage!='fixed':sc.strand('RTcDNA',[*rt_parts(),role('antisensetranscriptcDNA; actual firststrand')],label='originalRTprimer+newnormalDNAcDNA; possibleoriginal retained')
 if stage in ('dropletPCR','recovered','dnaindexed','rnaindexed'):
  d=[*cell_parts(),role('nativecopiedCS/gDNA targetchain; primer/endgeometryunassigned'),Segment('sourceR2Ntargetingoverhang',R2N)]
  r=[*cell_parts(),role('nativecopiedCSjunction; relativeplacementbounded'),role('retainedRTsamplebarcode',length_bp=8,molecular_feature=RT_SAMPLE),role('retainedN8UMI',length_bp=8,molecular_feature=UMI),role('targetedcDNA'),Segment('sourceRNAR2overhang',R2)]
  if stage in ('dnaindexed','rnaindexed'):
   selected=d if stage=='dnaindexed' else r
   selected.extend([role('installedi5sampleindex',length_bp=8,molecular_feature=I5),role('installedi7sampleindex',length_bp=8,molecular_feature=I7)])
  sc.strand('gDNAcopies',d,label='newnormalDNA gDNA-target class; '+('preferentiallyDNAindexamplified' if stage=='dnaindexed' else 'co-carriedtargetclass'))
  sc.strand('rnacopies',r,label='newnormalDNA cDNA-target class; '+('preferentiallyRNAindexamplified' if stage=='rnaindexed' else 'co-carriedtargetclass'))
 sc.strand('residualoligos',[role('independentfreeprimers/beadoligos/offtargetDNA orhydrolysates; noassumedphysicaldiscard')],label='independentpotentialcarryover')
 return sc

def oligos():
 return [f'<b>V2RT_primer1_8</b>:5′-{RNA_FORWARD}TCGCCTTANNNNNNNNTTTTTTTTTTTTTTTTVN-3′',f'<b>Droplet_amplification_fwd</b>:5′-{RNA_FORWARD}-3′',f'<b>gDNAforward</b>:5′-{CS}-gdna_specific-3′',f'<b>gDNAreverse</b>:5′-{R2N}-gdna_specific-3′',f'<b>Nextera_R1N_S502</b>:5′-{P5}CTCTCTAT{R1N_INDEX_ARM}-3′',f'<b>DNA_N701</b>:5′-{P7}TCGCCTTA{R2N_INDEX_ARM}-3′',f'<b>RNA_N701</b>:5′-{P7}TCGCCTTA{R2}-3′']
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def op(i,l,k,x,extra,state,reason=None):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state);g.react(i,l,k,((x,I.TEMPLATE),*extra),(y,),**({'reason':reason} if reason else {}),**({'changes_molecule':False} if k==K.GENERIC else {}));return y
 def clean(i,l,x,bounded=False):
  sel=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+'retainedfraction',M.MOLECULAR,state=x.state);w=g.product(i+'_waste','physicalnonretainedfraction; exactclassesunassigned',M.WASTE)
  g.react(i,l,K.MAGIC if bounded else K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w),**({'reason':'Sourceexplicitguidecleanup step; detailednativeconditions/partitionnotread in3-originalcap.'} if bounded else {}));g.discard(w,'Actualsourcephysicalcleanupfractions; modalitypurity notpresumed.');return y
 cell=g.input('fixed_cells','sourcefixed/permeabilizedcell input; selectedglyoxalroute',M.BIOLOGICAL,state=mixture('fixed'),primary=True)
 rtprimer=inp('rt_primer','sourceV2RT_primer1_8 representativefamily',M.OLIGO,ss('RTprimer',rt_parts()))
 rt=inp('rt_enzyme','MaximaHMinus;sourceprimaryEP0752; detailedmaterialsEP0751catalogdifference retainedunassigned',M.CATALYST)
 rtconditions=inp('rt_conditions','500µMdNTPeach;10min50°C,3cycles8/15/20/30/42/50°C,final50°C5min')
 x=op('rt','In situ barcodedpoly-dT reverse transcription',K.REVERSE_TRANSCRIPTION,cell,((rtprimer,I.STRUCTURAL),(rt,I.CATALYST),(rtconditions,I.CONDITION)),mixture('RT'))
 other=inp('other_rt_wells','independentRTwells pooledbeforeTapestri',M.MOLECULAR,mixture('RT'))
 pooled=g.product('pooled_rt','physicalpooledRTcells',M.MOLECULAR,state=x.state);g.react('pool_rt','PoolRTcells into10mlPBS/1%BSA',K.POOL,((x,I.TEMPLATE),(other,I.TEMPLATE)),(pooled,),changes_molecule=False)
 x=clean('rt_wash','500g5mincellpellet;resuspendMissionBioCellBuffer;count/dilute4000–8000cells/µl',pooled)
 device=inp('device','Tapestriv2 firstdropletencapsulation')
 enc=g.product('encapsulated','firstdropletpopulation; actualmixedcells',M.MOLECULAR,state=x.state)
 unused=g.product('unallocated_cells','possibleunallocatedsourcecells; amountunassigned',M.MOLECULAR,state=x.state)
 g.react('encapsulate','Encapsulate sourcecells in firstdropletpopulation',K.PARTITION,((x,I.TEMPLATE),(device,I.CONDITION)),(enc,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Unallocatedpopulation notassumedphysicallydestroyed.'));x=enc
 protease=inp('proteinase_k','sourceproteinaseK firstdropletlysis/digestion',M.CATALYST)
 x=op('lyse_digest','Lyse/digest firstdropletcells; nativeconditionsbounded',K.DEPROTEINIZATION,x,((protease,I.CATALYST),),x.state)
 f=inp('gdna_forward','CS–gdna_specific sourcepool',M.OLIGO,ss('gDNAforward',[Segment('sourceCS',CS),role('gDNAtargetspecific3arm')]))
 rf=inp('rna_forward','sourceDroplet_amplification_fwd',M.OLIGO,ss('rnaforward',[Segment('sourceV2forward',RNA_FORWARD)]))
 dr=inp('gdna_reverse','sourceR2N–gdna_specific',M.OLIGO,ss('gDNAreverse',[Segment('R2N',R2N),role('gDNAtargetspecific3arm')]))
 rr=inp('rna_reverse','sourceR2–cDNA_specific; panelmembersunassigned',M.OLIGO,ss('rnareverse',[Segment('knownRNAindexbindingprefix',R2),role('remainingnativeR2overhang/targetspecificarm; completepanelprimerunassigned')]))
 bead=inp('bead','MissionBiobarcodingbead; nativeorderedbases/photo-linkerunassigned',M.OLIGO,ss('beadnative',[role('sourcecellbarcode/nativeCScomplementcontext',molecular_feature=CELL1),role('secondcellbarcodecomponent',molecular_feature=CELL2)]))
 pcrpol=inp('droplet_polymerase','Tapestriv2 nativePCRpolymeraseunassigned',M.CATALYST)
 fuse=inp('fusion_conditions','Second dropletfusion addsforwardpools/PCRreagents/barcodingbead; nativefusionconditionsunassigned')
 x=op('fuse','Second droplet:addbarcodebead/forwardPCRmix toactualmixedlysate',K.GENERIC,x,((f,I.STRUCTURAL),(rf,I.STRUCTURAL),(dr,I.STRUCTURAL),(rr,I.STRUCTURAL),(bead,I.STRUCTURAL),(fuse,I.CONDITION)),x.state)
 uv=inp('uv_condition','sourceTapestriv2Step5UVtreatment; releaseidentity/conditionsunread')
 x=op('uv','UVtreat sourcebarcodingreaction; nativephotochemistrybounded',K.MAGIC,x,((uv,I.CONDITION),),x.state,reason='Directprotocol namesUVtreatment but nativeguide/linkerchemistryunread; no inventedphotolinker.')
 y=g.product('mixed_amplicons','jointcellbarcodedgDNA/cDNAamplicons pluspossibleoriginals',M.MOLECULAR,state=mixture('dropletPCR'))
 g.react('droplet_pcr','TargetedmultiplexPCR inactualmixedbarcodingdroplet',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(rf,I.STRUCTURAL),(dr,I.STRUCTURAL),(rr,I.STRUCTURAL),(bead,I.STRUCTURAL),(pcrpol,I.CATALYST)),(y,),geometry=PCRGeometry.unavailable(x.state,reason='ActualgDNA/ss-cDNA/RNA/nativebead mixture; panel-specific firstcycles/complete beadgeometryunassigned.'))
 breaker=inp('breaker','Sourceemulsionbreaking; nativereagent/orderdetailsunread')
 y=op('break_emulsion','Breaksourceemulsion beforePCRproductcleanup',K.MAGIC,y,((breaker,I.CONDITION),),mixture('recovered'),reason='Detailedprotocol explicitlybreaksemulsions; exactguide6nativeworkupunread.')
 y=clean('recover_amplicons','SourceTapestriv2Step6PCRproductcleanup',y,True)
 dna=g.product('dna_aliquot','15µl recoveredmixedamplicons forDNAindexPCR',M.MOLECULAR,state=y.state)
 rna=g.product('rna_aliquot','15µl recoveredmixedamplicons forRNAindexPCR',M.MOLECULAR,state=y.state)
 rem=g.product('remaining_amplicons','unallocatedrecoveredmixedamplicons',M.MOLECULAR,state=y.state)
 rule=inp('split_rule','15µl inputeachindexPCR; allocationdoesnotpurifymodalities')
 g.react('split_index_inputs','Takephysicalmixedaliquots fordistinctoverhangPCRs',K.SAMPLING,((y,I.TEMPLATE),(rule,I.CONDITION)),(dna,rna,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Unallocatedmixedmaterialnotassumeddestroyed.'))
 i5=inp('index_i5','V2Nextera_R1N_S502 representative8ntfamily',M.OLIGO,ss('i5',[Segment('P5',P5),Segment('sourcei5example','CTCTCTAT',feature=I5),Segment('R1Narm',R1N_INDEX_ARM)]))
 di7=inp('dna_i7','V2Nextera_R2N_N701 representative8ntfamily',M.OLIGO,ss('dnai7',[Segment('P7',P7),Segment('sourcei7example','TCGCCTTA',feature=I7),Segment('R2Narm',R2N_INDEX_ARM)]))
 ri7=inp('rna_i7','V2Illumina_R2_N701; selectR2 pertext/Table2, notdetailedRNA tableR2N',M.OLIGO,ss('rnai7',[Segment('P7',P7),Segment('sourcei7example','TCGCCTTA',feature=I7),Segment('R2arm',R2)]))
 mix=inp('library_polymerase','MissionBioLibraryPCRMM2×;source10cycles98/62/72°C',M.CATALYST)
 for name,input_,reverse,stage,lens in [('dna',dna,di7,'dnaindexed',{'Read 1':161,'Read 2':161,'Index 1 (i7)':8,'Index 2 (i5)':8}),('rna',rna,ri7,'rnaindexed',{'Read 1':76,'Read 2':46,'Index 1 (i7)':8,'Index 2 (i5)':8})]:
  output=g.product(name+'_indexed',name+'preferentialindexedcopies withpossibleotherclasses',M.MOLECULAR,state=mixture(stage))
  g.react(name+'_index_pcr',name+'overhang-specificPCR onactualmixed15µl input;10cycles',K.PCR,((input_,I.TEMPLATE),(i5,I.STRUCTURAL),(reverse,I.STRUCTURAL),(mix,I.CATALYST)),(output,),geometry=PCRGeometry.unavailable(input_.state,reason='Actualrecoveredmixture; primerchoice supplies preferentialamplification, no modalitypurificationclaim.'))
  output=clean(name+'_cleanup','Sourcev2Guide7.7librarycleanup/8quantify-normalize',output,True)
  lib=g.product(name+'_library',name+'sequencinglibrary; co-recoveredclassesnotassumedabsent',M.LIBRARY,state=output.state)
  g.react(name+'_retain','Retain'+name+'nativeindexedlibrary',K.GENERIC,((output,I.TEMPLATE),),(lib,),changes_molecule=False)
  g.sequence(lib,SequencingPlan(read_lengths=lens,acquired_roles=tuple(lens),limitation=SEQUENCING_ENDING+' Nativeprimer/targetgeometrybounded.'))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('fixed_cells','rt_primer','rt_enzyme','rt_conditions','rt','other_rt_wells','pool_rt','rt_wash_selection','rt_wash','device','encapsulate','proteinase_k','lyse_digest','gdna_forward','rna_forward','gdna_reverse','rna_reverse','bead','droplet_polymerase','fusion_conditions','fuse','uv_condition','uv','droplet_pcr','breaker','break_emulsion','recover_amplicons_selection','recover_amplicons','split_rule','split_index_inputs','index_i5','dna_i7','rna_i7','library_polymerase','dna_index_pcr','rna_index_pcr','dna_cleanup_selection','dna_cleanup','rna_cleanup_selection','rna_cleanup','dna_retain','rna_retain','dna_library','rna_library')
