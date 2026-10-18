"""LiMCA2024 native compartment-specific chemistry, cited pex2423v1 scope explicit."""
from chemdraw import Construct,Scene,Segment,feature,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,Terminal,Disposition,SequencingPlan
from pcr import PCRGeometry
import seqprimers as sp
TITLE='LiMCA — source-separated cytoplasm and nuclear contacts'
NOTES='01_limca.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-024-02239-0">Wu et al., Nature Methods2024</a>; directly cited pex2423v1 protocol and its original attachedLiMCAOligos workbook.'
SUMMARY='Physical5µl cytoplasmic supernatant and nuclear pellet from one cell follow distinct Smart-seq2/TruePrep and NlaIII/META20 routes. Matching uses experimental tube labels, not a shared molecular cell barcode.'
CAVEAT='Version scope: primary2024Methods oligo-dT preheat72°C5min and RNA0.6×/0.2×selection are selected; citedauthorprotocolv1 says70°C5min and0.55×/0.25×. These ordinary source-version differences are not silently merged. The v1 workbook supplies native Smartseq2/META/UDP components. Printed META20transfer strands end in13ntME sequence while the companion ME is19nt; the drawing preserves13bppairing and6ntcompanion overhang, without appending missing bases or proving loaded transposase geometry. Native META insertion/PCRfirstcycle/DM5xx-DM7xx indexassembly and selectedcellmembers remain bounded. Molecularindex widths are separate from acquisition; no generic8nt shell.'
RT_HANDLE='AAGCAGTGGTATCAACGCAGAGTAC'
TSO_DNA='AAGCAGTGGTATCAACGCAGAGTACAT'
ISPCR='AAGCAGTGGTATCAACGCAGAGT'
META20_1='AGAAGCCGTGTGCCGGTCTAAGATGTGTATAAG'
ME19='CTGTCTCTTATACACATCT'
META40_1='ACACTCTTTCCCTACACGACGCTCTTCCGATCT'+META20_1
META40_2='GACTGGAGTTCAGACGTGTGCTCTTCCGATCT'+META20_1
UDP5='AATGATACGGCGACCACCGAGATCTACACTCGTGGAGCGTCGTCGGCAGCGTC'
UDP7='CAAGCAGAAGACGGCATACGAGATCGCTCAGTTCGTCTCGTGGGCTCGG'
META_FEATURE=feature('limca.meta_family','inline_barcode','whitelist',whitelist='source20META amplification tagfamily',note='technical amplification tag family; not a sourcecellbarcode or UMI claim')
CELL_INDEX=feature('limca.cell_library','sample_index','unknown',note='percell tube/library indexmember, DNAindexwidth unassigned')
RNA5=feature('limca.rna_i5','sample_index','whitelist',whitelist='sourceUDP family',note='UDP0001example10ntphysicalcode; acquisitionunassigned')
RNA7=feature('limca.rna_i7','sample_index','whitelist',whitelist='sourceUDP family',note='UDP0001example10ntphysicalcode; acquisitionunassigned')
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_ENDING='Published LiMCA PE150 onHiSeqX10/NovaSeq; native index acquisition unassigned. DNA/RNA tube labels link the cell, not molecular concatenation.'
READ_SPECS=tuple(sp.read_spec(r,150 if r.startswith('Read') else 'unassigned','Contact/transcript insert' if r.startswith('Read') else 'Native index acquisition allocation unassigned','include' if r.startswith('Read') else 'undetermined') for r in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
def rt_parts():return [Segment('source RT handle',RT_HANDLE),Segment('dT30','T'*30),Segment('VNanchor','VN',placeholder=True)]
def ss_cdna(switched=False):
 sc=Scene();parts=[*rt_parts(),role('transcript antisense cDNA; native extent unassigned')]
 if switched:parts.extend([role('native switch junction/terminalGcopy extent unassigned'),Segment('copied TSO DNAprefix',revcomp(TSO_DNA))])
 sc.strand('cDNA',parts,label='actual '+('switched' if switched else 'first')+' ssDNA precursor')
 sc.strand('originals',[role('possible originalRNA/cytoplasmicDNA/oligos; recovery unresolved')],label='independent residualorigins');return sc

def native_state(kind,stage,indexed=False):
 sc=Scene()
 parts=[role('native '+stage+' endA'),role(kind+' DNA/contact/copy extent unassigned')]
 if kind=='nuclear':parts.insert(1,role('sourceMETA tagfamily; member/orientation unassigned',molecular_feature=META_FEATURE))
 if indexed:
  parts.extend([role('native i5libraryindex',molecular_feature=RNA5 if kind=='RNA' else CELL_INDEX),role('native i7libraryindex',molecular_feature=RNA7 if kind=='RNA' else CELL_INDEX)])
 sc.strand('native DNA',parts,label='actual '+kind+' population; '+stage)
 sc.strand('residuals',[role('possible sourceoriginalRNA/DNA, native primers, priorcopies and offproducts; physical recovery unassigned')],label='independent possible residual lineages');return sc

def cell_scene():
 sc=Scene();sc.strand('RNA',[role('cytoplasmicRNA population')],label='RNA in cell');sc.strand('DNA',[role('nuclear genomicDNA')],label='independent nuclear DNA');return sc

def fraction_scene(cytoplasm):
 sc=Scene();sc.strand('primary',[role('cytoplasmicRNA' if cytoplasm else 'nuclear genomicDNA')],label='source-enriched physical compartment')
 sc.strand('possible carryover',[role('possible cytoplasmicDNA/other material' if cytoplasm else 'possible retained nuclearRNA/other material')],label='independent residual population; amount/recovery unassigned')
 return sc

def cut_scene():
 sc=Scene()
 for k in ('A','B'):
  body=role(k+' NlaIIIfragment');tail=Segment(k+'3′CATG','CATG')
  sc.strand(k,[body,tail],mod5='p',label='NlaIII3′CATG end')
  opp=complement_segments([body]);sc.anneal(k+' opposite',[*opp,Segment(k+'opposed3′CATG','CATG')],to=k,pair=(opp[0].name,body.name),unpaired=(opp[0].name,),mod5='p')
  sc.overhang(k,tail.name,polarity='3-prime',owner='substrate',text='sourceNlaIII3′CATG')
 sc.strand('possible RNA',[role('possible retained nuclearRNA/other original material; recovery unassigned')],label='independent residualorigins')
 return sc

def contact_scene():
 parts=[role('locusA'),Segment('NlaIIIcohesive contact junction','CATG'),role('locusB')]
 sc=Scene.duplex(parts,label='original nuclear3CcontactDNA; no biotin fill',unpaired=tuple(p.name+"'" for p in parts if p.is_role_token()))
 sc.strand('possible RNA',[role('possible retained originalRNA/protein/other material; recovery unresolved')],label='independent residualorigins');return sc

def meta_adapter():
 core=META20_1[-13:]
 if revcomp(ME19[6:])!=core:raise ValueError('sourceMETA13bp overlap mismatch')
 sc=Scene();sc.strand('transfer',[Segment('METAamplifier tag',META20_1[:-13],feature=META_FEATURE),Segment('printedME13',core)],label='actual printedMETA20_1')
 sc.anneal('ME19',[Segment('ME5′overhang6',ME19[:6]),Segment('pairedME13',ME19[6:])],to='transfer',pair=('pairedME13','printedME13'),unpaired=('ME5′overhang6',),label='printed19bpME; no silent topstrandcompletion')
 return sc
RNA_LIB=Construct([role('native RNAendA'),role('transcriptinsert'),role('native UDPi5',length_bp=10,molecular_feature=RNA5),role('native UDPi7',length_bp=10,molecular_feature=RNA7)],name='LiMCA RNA native endpoint')
DNA_LIB=Construct([role('native META/DMendA'),role('META amplificationtag',molecular_feature=META_FEATURE),role('NlaIIIcontactinsert'),role('native cell-libraryindexregion',molecular_feature=CELL_INDEX)],name='LiMCA nuclear native endpoint')
FINAL_LIBRARY=None
FINAL_LIBRARIES=(('RNA native reference',RNA_LIB,(),'Source components; complete nativeassembly unassigned.',SEQUENCING_ENDING),('Nuclear native reference',DNA_LIB,(),'Source META/DM assembly unassigned.',SEQUENCING_ENDING))
def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def op(i,l,k,x,extra,state):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state);g.react(i,l,k,((x,I.TEMPLATE),*extra),(y,));return y
 def clean(i,l,x):
  sel=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+'recovered',x.kind,state=x.state);w=g.product(i+'_waste','physical nonretainedfraction',M.WASTE)
  g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));g.discard(w,'Source physical fraction; class purity/recovery not presumed.');return y
 def pcr(i,l,x,f,r,pol,state,symmetric=False):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state);g.react(i,l,K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),*(((r,I.STRUCTURAL),) if r else ()),(pol,I.CATALYST)),(y,),geometry=PCRGeometry.unavailable(x.state,reason='Actual incoming native '+i+' precursor; full strand/primer geometry bounded, not finalendpoint substituted.',symmetric_species=symmetric));return y
 cell=g.input('cell','one picked viable cell',M.BIOLOGICAL,state=cell_scene(),primary=True)
 soft=inp('soft_lysis','7µl0.45%IGEPAL30minice;vortex;500g5min4°C')
 cyt=g.product('cytoplasm','5µlsupernatant fraction; RNAenriched not provenpure',M.MOLECULAR,state=fraction_scene(True));nuc=g.product('nucleus','remaining2µlnuclearpellet; possibleRNAcarryover',M.MOLECULAR,state=fraction_scene(False))
 g.react('separate','Physicallytransfer5µlsupernatant aftersoftlysis/spin',K.SPLIT,((cell,I.TEMPLATE),(soft,I.CONDITION)),(cyt,nuc),changes_molecule=False)
 dt=inp('rt_primer','sourceoligo-dT30VN; preheat72°CprimaryMethods',M.OLIGO,Construct(rt_parts()));rt=inp('ssii','SSII/betaine/Mg/GTP;42°C90min+10×50/42°C',M.CATALYST)
 first=op('rt','Reverse-transcribe cytoplasmicRNA',K.REVERSE_TRANSCRIPTION,cyt,((dt,I.STRUCTURAL),(rt,I.CATALYST)),ss_cdna())
 tso_scene=Scene();tso_scene.strand('TSO',[Segment('TSODNAprefix',TSO_DNA),Segment('RNA/LNA GGG','GGG',note='source rGrG+G')]);tso=inp('tso','sourceTSO-LNA',M.OLIGO,tso_scene)
 switched=op('switch','Template-switch in sameRTreaction',K.TEMPLATE_SWITCHING,first,((tso,I.STRUCTURAL),(rt,I.CATALYST)),ss_cdna(True))
 isp=inp('ispcr','sourceISPCRsingle species',M.OLIGO,ISPCR);kapa=inp('kapa','KAPAHiFiHotStart21cycles',M.CATALYST)
 rna=pcr('rna_preamp','Amplify switchedssDNA21cycles',switched,isp,None,kapa,native_state('RNA','preamp'),True);rna=clean('rna_pre_cleanup','0.7×AMPure cleanup',rna)
 tag=inp('rna_loaded','nativeVazymeTD502loadedends; strandsunassigned',M.ADAPTER);tn=inp('rna_tn5','TruePrepTTE/TTBL55°C10min;selected5ng',M.CATALYST)
 rna=op('rna_tag','Tagment actual cDNAamplicons',K.TAGMENTATION,rna,((tag,I.STRUCTURAL),(tn,I.CATALYST)),native_state('RNA','tagged'))
 stop=inp('rna_stop','2.5µlTSstopbuffer');rna=op('rna_stop_step','Stop source transposition',K.QUENCHING,rna,((stop,I.CONDITION),),rna.state)
 f=inp('rna_i5','UDP0001i5 familyexample; selectedmemberunassigned',M.OLIGO,UDP5);r=inp('rna_i7','UDP0001i7 familyexample; selectedmemberunassigned',M.OLIGO,UDP7);q=inp('rna_q5','sourceQ5mastermix',M.CATALYST)
 rna=op('rna_extend','Authorv1native72°C5min initialextension',K.EXTENSION,rna,((f,I.STRUCTURAL),(r,I.STRUCTURAL),(q,I.CATALYST)),native_state('RNA','extended'))
 rna=pcr('rna_index','Authorv1RNAindexPCR9cycles',rna,f,r,q,native_state('RNA','indexed',True))
 others=inp('rna_others','other independentlyindexedRNAcelllibraries',M.MOLECULAR,native_state('RNA','indexed',True));pool=g.product('rna_pool','pooledRNAcelllibraries',M.MOLECULAR,state=rna.state);g.react('rna_pooling','PoolRNAcelllibraries',K.POOL,((rna,I.TEMPLATE),(others,I.TEMPLATE)),(pool,),changes_molecule=False)
 rna=clean('rna_dcc','SourceDCC5purification',pool);rna=clean('rna_size1','Primary2024RNA0.6×selection',rna);rna=clean('rna_size2','Primary2024RNA0.2×recovery',rna)
 rf=g.product('rna_library',RNA_LIB.name,M.LIBRARY,state=rna.state);g.react('rna_retain','RetainRNAfraction',K.GENERIC,((rna,I.TEMPLATE),),(rf,),changes_molecule=False)
 # Nuclear3C then source-defined META20/META40/index operations.
 fix=inp('pfa','8µl2.5%PFA10minRT;glycinequench');nuc=op('fix','Fix isolatednucleus',K.CROSS_LINKING,nuc,((fix,I.CONDITION),),nuc.state);nuc=clean('nuc_wash','Pellet/wash/remove17µlsupernatants',nuc)
 perm=inp('nuc_permeabilization','0.3%SDS62°C10min;Triton37°C15min');nuc=op('permeabilize','Permeabilize nuclearpellet',K.PERMEABILIZATION,nuc,((perm,I.CONDITION),),nuc.state)
 nl=inp('nlaiii','NlaIII37°C2h',M.CATALYST);nuc=op('digest','DigestNlaIII3′CATGends',K.DIGESTION,nuc,((nl,I.CATALYST),),cut_scene());nuc=clean('digest_wash','Pellet/remove17µl;T4bufferwash',nuc)
 lig=inp('ligase','T4DNA ligaseRT2.5h',M.CATALYST);nuc=op('ligate','Ligate original3Ccontacts without fill-in',K.PROXIMITY_LIGATION,nuc,((lig,I.CATALYST),),contact_scene());nuc=clean('ligate_wash','Pellet/remove18µlsupernatant',nuc)
 prot=inp('qp','QIAGENprotease/DTTTriton50°C1h65°C1h70°C15min',M.CATALYST);nuc=op('lyse','Lyse nuclearcontact sample',K.DEPROTEINIZATION,nuc,((prot,I.CATALYST),),contact_scene())
 meta=inp('meta_adapter','equal20sourceMETAtransferoligos annealedwithME19; example13bpoverlap',M.ADAPTER,meta_adapter());tn=inp('meta_tn5','VazymeTn5assembledRT30min;transposition55°C10min',M.CATALYST)
 nuc=op('meta_tag','Tagment source3Clysate withMETA20',K.TAGMENTATION,nuc,((meta,I.STRUCTURAL),(tn,I.CATALYST)),native_state('nuclear','META tagged'))
 remove=inp('meta_protease','NaCl/EDTA/QIAGENprotease50°C30min70°C15min',M.CATALYST);nuc=op('meta_remove','RemoveTn5protein byprotease',K.DEPROTEINIZATION,nuc,((remove,I.CATALYST),),nuc.state)
 m20=inp('meta20_primers','sourceMETA20 primerpool; exampleMETA20-1',M.OLIGO,META20_1);q=inp('meta_q5','Q5/MgCl2amplificationmix',M.CATALYST)
 nuc=op('meta_extend','72°C5min precycleextension',K.EXTENSION,nuc,((m20,I.STRUCTURAL),(q,I.CATALYST)),native_state('nuclear','META extended'))
 nuc=pcr('meta_pcr','META20preamp13cycles',nuc,m20,None,q,native_state('nuclear','META preamp'),True)
 for n in (1,2):
  ex=inp('exo'+str(n),'ExoI37°C30min80°C20min',M.CATALYST);nuc=op('exo_step'+str(n),'Digest excessprimers; productsnotphysicaldiscard',K.DIGESTION,nuc,((ex,I.CATALYST),),nuc.state)
  if n==1:
   mf=inp('meta40_f','META40-1example sourcepool',M.OLIGO,META40_1);mr=inp('meta40_r','META40-2example sourcepool',M.OLIGO,META40_2);nuc=pcr('meta40_pcr','IncorporateMETA40adapterarms2cycles',nuc,mf,mr,q,native_state('nuclear','META40'))
  else:
   df=inp('dm5','DM5xxnativeindexprimer; fullsource/memberunassigned',M.OLIGO);dr=inp('dm7','DM7xxnativeindexprimer; fullsource/memberunassigned',M.OLIGO);nuc=pcr('dna_index','DM5/DM7libraryPCR2cycles',nuc,df,dr,q,native_state('nuclear','indexed',True))
 nuc=clean('dna_dcc','PerwellDCC5elute20µl',nuc);others=inp('dna_others','other independentlyindexednuclear libraries',M.MOLECULAR,native_state('nuclear','indexed',True));pool=g.product('dna_pool','pooled nuclear libraries',M.MOLECULAR,state=nuc.state);g.react('dna_pooling','Pool independently indexednuclear libraries',K.POOL,((nuc,I.TEMPLATE),(others,I.TEMPLATE)),(pool,),changes_molecule=False);nuc=clean('dna_cleanup','0.8×AMPure afterpooling',pool)
 df=g.product('dna_library',DNA_LIB.name,M.LIBRARY,state=nuc.state);g.react('dna_retain','Retainnuclearfraction',K.GENERIC,((nuc,I.TEMPLATE),),(df,),changes_molecule=False)
 for x,lib in ((rf,RNA_LIB),(df,DNA_LIB)):g.sequence(x,SequencingPlan(library=lib,acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('cell','soft_lysis','separate','rt_primer','ssii','rt','tso','switch','ispcr','kapa','rna_preamp','rna_pre_cleanup','rna_pre_cleanup_selection','rna_loaded','rna_tn5','rna_tag','rna_stop','rna_stop_step','rna_i5','rna_i7','rna_q5','rna_extend','rna_index','rna_others','rna_pooling','rna_dcc','rna_dcc_selection','rna_size1','rna_size1_selection','rna_size2','rna_size2_selection','rna_retain','rna_library','pfa','fix','nuc_wash','nuc_wash_selection','nuc_permeabilization','permeabilize','nlaiii','digest','digest_wash','digest_wash_selection','ligase','ligate','ligate_wash','ligate_wash_selection','qp','lyse','meta_adapter','meta_tn5','meta_tag','meta_protease','meta_remove','meta20_primers','meta_q5','meta_extend','meta_pcr','exo1','exo_step1','meta40_f','meta40_r','meta40_pcr','exo2','exo_step2','dm5','dm7','dna_index','dna_dcc','dna_dcc_selection','dna_others','dna_pooling','dna_cleanup','dna_cleanup_selection','dna_retain','dna_library')
def oligos():return ('oligo-dT—5′-'+RT_HANDLE+'T30VN-3′','TSO-LNA—5′-'+TSO_DNA+'rGrG+G-3′','ISPCR—5′-'+ISPCR+'-3′','META20-1—5′-'+META20_1+'-3′','19bpME—5′-'+ME19+'-3′','META40-1—5′-'+META40_1+'-3′','META40-2—5′-'+META40_2+'-3′','UDP0001i5—5′-'+UDP5+'-3′','UDP0001i7—5′-'+UDP7+'-3′')
