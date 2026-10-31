"""Published SIMPLE-seq physical preparation with locally scoped molecular interfaces.

Occupied genomic strands and study-wide barcode inventories are unassembled. No
unread predecessor recipe or authors' later bench variant supplies hidden defaults.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import illumina,nextera
AUTHORITY='Bai et al., published Nature Biotechnology DOI10.1038/s41587-024-02148-9, Methods2582–2663/2671–2820; original SupplementaryTable1_1 rows3–12. Main governs; authors bench/GEO alternatives not merged.'
LINKER1='CTCGAATGCCTCAAGCACGTGGAT'
LINKER2='GGTCTGAGTTCGCACCGAAACATCGGCCAC'
BLOCK1='ATCCACGTGCTTGAGAGGCCAGAGCATTCG'
BLOCK2='GTGGCCGATGTTTCGGTGCGAACTCAGACC'
ME=nextera.ME
BARCODED='GCATTCGAGACGCA'+ME
BRIDGE_TOP=revcomp(nextera.ADAPTOR_S5+'CTGACG')
BRIDGE_FIXED=ME+'CTGACG'
INDICATOR=nextera.ADAPTOR_S5+'CTGACG'
UNIVERSAL=illumina.P5+nextera.S5
INDEX_LEFT=illumina.P7
INDEX_RIGHT=illumina.TRUSEQ_READ2[:-1]
RAW_ORDERINGS={
 'Linker-R01':LINKER1,'Linker-R02':LINKER2,'Blocker-R01':BLOCK1,'Blocker-R02':BLOCK2,
 'pMENTs':'5Phos/CTGTCTCTTATACACATCT/ddC/',
 'Bridge Adaptor-bottom':'AGATGTGTATAAGAGACAGCTGACGNNNNNN/-NH2',
 'Bridge Adaptor-top':'5Phos/CGTCAGCTGTCTCTTATACACATCTGACGCTGCCGACGA',
 'P7-indicator primer':'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG/caC/TGA/caC/G',
 'P7 index primer':'CAAGCAGAAGACGGCATACGAGAT[6nt i7 index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC',
 'N5 Universal primer':UNIVERSAL,
}
def single(seq,label,mod5='',mod3=''):
 sc=Scene();sc.strand(label,[Segment(label,seq)],mod5=mod5,mod3=mod3);return sc

def transposon():
 sc=Scene();sc.strand('barcoded DNA',[Segment('ligation overhang','GCATTCGAG'),Segment('sample barcode','ACGCA',feature=feature('simple.tn5','cell_barcode','fixed',group='simple.cell',part='tn5_component',note='Representative source DNA01; one component of joint cell identity, not unique per-cell or sample i7')),Segment('ME',ME)],mod5='P')
 sc.anneal('pMENTs',[Segment('ME reverse',revcomp(ME)),Segment('additional terminal ddC','C')],to='barcoded DNA',pair=('ME reverse','ME'),unpaired=('additional terminal ddC',),mod5='P',mod3='ddC')
 return sc

def bridge():
 sc=Scene();sc.strand('bridge top',[Segment('paired core',BRIDGE_TOP[:25]),Segment('PCR tail',BRIDGE_TOP[25:])],mod5='P')
 sc.anneal('bridge bottom',[Segment('paired core reverse',BRIDGE_FIXED),Segment('random capture overhang','NNNNNN')],to='bridge top',pair=('paired core reverse','paired core'),mod3='NH2')
 sc.ligation_end('bridge top','paired core',kind='splinted',polarity='5-prime',owner='adapter')
 return sc

def indicator():
 sc=Scene();sc.strand('indicator',[Segment('common handle',INDICATOR[:-6]),Segment('indicator1','C'),Segment('intervening','TGA'),Segment('indicator2','C'),Segment('terminal','G')]);return sc

def index_primer():
 sc=Scene();sc.strand('P7 index primer',[Segment('P7',INDEX_LEFT),Segment('source6nt index','[i7]',placeholder=True,length_bp=6,feature=feature('simple.i7','sample_index','unknown',note='Source6nt member value unassigned')),Segment('primer right',INDEX_RIGHT)]);return sc

class B:
 def __init__(self):
  self.g=ProtocolGraph('SIMPLE-seq — published nuclei preparation',authority=AUTHORITY,representation='published_instance');self.x=self.g.input('cells','Source cells, with native genomic5mC/5hmC; no chosen genomic sequence',K.BIOLOGICAL,primary=True)
 def stock(self,k,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(k,label,kind,**kw),role)
 def cond(self,k,label):return self.stock(k,label,K.CONDITION,R.CONDITION)
 def op(self,k,label,kind=S.MAGIC,ins=(),product=None,chemical=None,**kw):
  y=self.g.product(k+'_product',product or 'Actual source preparation; whole occupied strands, recovery, specificity and conversion fractions unassembled',K.MOLECULAR,chemical_state=chemical)
  if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; full occupied molecules/intermediates and quantitative products unassembled, not absence evidence.')
  if kind==S.GENERIC:kw['changes_molecule']=False
  if kind in {S.CROSS_LINKING,S.CHEMICAL_CONVERSION} and not any(r==R.CONDITION for _,r in ins):
   ins=tuple(ins)+(self.cond(k+'_conditions',label),)
  self.g.react(k,label,kind,((self.x,R.TEMPLATE),)+tuple(ins),(y,),**kw);self.x=y;return y
 def clean(self,k,label):
  y=self.g.product(k+'_kept','Retained source fraction; molecular purity/carryover unassigned',K.MOLECULAR);w=self.g.product(k+'_removed','Actual excluded workup fraction',K.WASTE)
  self.g.react(k,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(k+'_conditions',label)),(y,w));self.g.discard(w,'Physical removed source workup, no ideal separation.');self.x=y
 def select(self,k,label):
  y=self.g.product(k+'_selected',label,K.MOLECULAR);rest=self.g.product(k+'_remainder','Actual untransferred source stock; later fate unassigned',K.MOLECULAR)
  self.g.react(k,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,rest));self.g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred source stock retained.'));self.x=y

def barcode_supplies(g,key,label):
 stock=g.input(key+'_supplies',label+'; aggregate actual supplies, stockcount/packaging unassigned',K.MOLECULAR)
 ol=g.product(key+'_oligos','Copresent source barcode classes; whole member bodies unexpanded',K.OLIGO)
 linker=g.product(key+'_linker','Copresent actual source linker fromsame supplies',K.OLIGO,state=single(LINKER1 if key=='ligation1' else LINKER2,'SourceLinkerR01' if key=='ligation1' else 'SourceLinkerR02'))
 lig=g.product(key+'_ligase','Actual source T4DNA ligase component',K.CATALYST)
 bsa=g.product(key+'_bsa','Actual source BSA molecular component',K.MOLECULAR)
 buf=g.product(key+'_buffer','Actual source ligation buffer component',K.CONDITION)
 g.react(key+'_components','Actual supply components, notextra fresh reagent additions',S.GENERIC,((stock,R.TEMPLATE),),(ol,linker,lig,bsa,buf),changes_molecule=False)
 return((ol,R.STRUCTURAL),(linker,R.STRUCTURAL),(lig,R.CATALYST),(bsa,R.STRUCTURAL),(buf,R.CONDITION))

# Fixed reviewed declarations follow the selected published Methods and actual workbook,
# independently of graph traversal. This is preparation scope, not an absence audit.
BIO_SCOPE=('cells','fixative','fix_conditions','fix','glycine','quench','cell_wash_conditions','cell_wash','nib','nib_conditions','isolate','nuclei_wash_conditions','nuclei_wash','sds','sds_conditions','deplete','triton','triton_conditions','quench_sds','nuclei_pellet_conditions','nuclei_pellet')
BARCODE_SCOPE=('tn5_oligos','tn5_oligo_components','tn5_top','tn5_bottom','tn5_anneal_conditions','tn5_anneal','tn5_enzyme','tn5_load_conditions','tn5_load','tn5_allocate','other_tn5','tag_components','tag_enzyme','tag_representative_adapter','tag_other_adapters','tag_other_material','tag_conditions','tag_array','edta','stop_tag','pool_tag_conditions','pool_tag','ligation1_supplies','ligation1_components','ligation1_oligos','ligation1_linker','ligation1_ligase','ligation1_bsa','ligation1_buffer','ligation1_conditions','barcode1','block1','block1_components','block1_oligo','block1_buffer','block1_conditions','block_round1','pool1_conditions','pool1','ligation2_supplies','ligation2_components','ligation2_oligos','ligation2_linker','ligation2_ligase','ligation2_bsa','ligation2_buffer','ligation2_conditions','barcode2','termination','termination_components','termination_oligo','termination_edta','stop_barcode','nuclei_count_conditions','nuclei_count','sublibrary','lysis_supplies','lysis_protease','lysis_conditions','lyse','dna_clean_conditions','dna_clean')
FINISH_SCOPE=('denature1_conditions','denature1','exoi','exoi_conditions','exoi1','exo_clean_conditions','exo_clean','denature2_conditions','denature2','bridge_adapter','bridge_ligase','bridge_peg','bridge_conditions','bridge_ligate','bridge_clean_conditions','bridge_clean','spike','add_spike','preoxidation_beads_conditions','preoxidation_beads','p6_conditions','p6','naoh','denature_chemical_conditions','denature_chemical','ruthenate','oxidation_conditions','oxidize_hmc','oxide_p6_conditions','oxide_p6','malononitrile','malono_conditions','label_hmc','indicator_primer','extension_nt','robust','extension_conditions','record_hmc','exoi2_enzyme','exoi2_conditions','exoi2','extension_clean_conditions','extension_clean','tet_enzyme','tet_cosubstrates','tet_conditions','oxidize_mc','protease','protease_conditions','proteolyse','tet_clean_conditions','tet_clean','borane','borane_conditions','reduce_mc','borane_clean_conditions','borane_clean','pcr_supplies','pcr_components','pcr_forward','pcr_reverse','pcr_enzyme','pcr_nt','pcr_other','pcr_conditions','index_pcr','final_clean_conditions','final_clean','loading','readiness','library')

def protocol_graph():
 b=B();g=b.g
 b.op('fix','Published1percent formaldehyde10minRT',S.CROSS_LINKING,(b.stock('fixative','Actual406ul37percentformaldehyde to10mlcell suspension'),))
 b.op('quench','Published800ul2.5Mglycine then5minice',S.GENERIC,(b.stock('glycine','Actual glycine quench reagent'),))
 b.clean('cell_wash','550g8min4C pellet/PBS wash/pellet')
 b.op('isolate','Source NIB20minice nuclei isolation',S.EXTRACTION,(b.stock('nib','Actual NIB inventory20mMHEPES/10mMNaCl/3mMMgCl2/.1percentIgepal/protease inhibitors'),b.cond('nib_conditions','Source NIB20minice')))
 b.clean('nuclei_wash','500g5min4C nuclei recovery')
 b.op('deplete','Source.3percentSDS60C10min',S.PERMEABILIZATION,(b.stock('sds','Actual SDS inNEBuffer2.1'),b.cond('sds_conditions','Source60C850rpm10min')))
 b.op('quench_sds','Source TritonX10037C60min',S.GENERIC,(b.stock('triton','Actual200ul10percentTritonX100'),b.cond('triton_conditions','Source37C850rpm60min')))
 b.clean('nuclei_pellet','Source100g4C10min nucleus pellet')
 # A representative source ordering; other eleven barcode stocks remain copresent,
 # not silently replaced by this one example.
 ol=g.input('tn5_oligos','Aggregateactual source DNA01 andpMENTs oligo contributions; separate stockpackaging/volumesunassigned',K.MOLECULAR)
 top=g.product('tn5_top','Actual source5PbarcodedDNA01 strand',K.OLIGO,state=single(BARCODED,'DNA01',mod5='P'))
 bottom=g.product('tn5_bottom','Actual source5P MEcomplement plusadditionalddC',K.OLIGO,state=single(revcomp(ME)+'C','pMENTs',mod5='P',mod3='ddC'),chemical_state={'additional_terminal_modified_base':(len(ME),'ddC'),'source_ordering':RAW_ORDERINGS['pMENTs']})
 g.react('tn5_oligo_components','Actualoligo inventory bookkeeping, notanotherfreshaddition',S.GENERIC,((ol,R.TEMPLATE),),(top,bottom),changes_molecule=False)
 duplex=g.product('tn5_duplex','Same annealed source duplex; complete stock stoichiometry unassigned',K.ADAPTER,state=transposon())
 g.react('tn5_anneal','Source95C5min,1C/mincool to4C',S.HYBRIDIZATION,((top,R.TEMPLATE),(bottom,R.STRUCTURAL),b.cond('tn5_anneal_conditions','Published source temperature ramp')),(duplex,))
 loaded=g.product('loaded_tn5','Actual representative loaded-transposome stock, source9ul10uMduplex+10ul.5mgmlTn5+17ulcouplingbuffer',K.MOLECULAR)
 g.react('tn5_load','Source30C30min Tn5 loading',S.MAGIC,((duplex,R.TEMPLATE),b.stock('tn5_enzyme','Actual VazymeS601 unloaded Tn5',K.CATALYST,R.CATALYST),b.cond('tn5_load_conditions','Source30C30min; actual source coupling buffer')),(loaded,),reason='Loaded complex stoichiometry/free companions unassembled; no invented exact occupied transposome.')
 selected=g.product('tn5_selected','Actual representative loaded-stock allocation; volume unassigned',K.MOLECULAR);rest=g.product('tn5_remainder','Untransferred loaded stock; exact amount unassigned',K.MOLECULAR)
 g.react('tn5_allocate','Select actual loaded-stock contribution',S.SPLIT,((loaded,R.TEMPLATE),),(selected,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred source factory stock.'))
 other=g.input('other_tn5','Actual source other11barcoded loaded-transposome stockclasses; original whole member orderings/volumes unexpanded',K.MOLECULAR)
 enzyme=g.product('tag_enzyme','Copresent Tn5 catalyst components ofsame source loaded-stock contributions',K.CATALYST)
 representative=g.product('tag_representative_adapter','Copresent representative DNA01/pMENTs component ofselected loadedstock, notfresh freeadapter',K.ADAPTER,state=transposon())
 adapters=g.product('tag_other_adapters','Copresent other11barcode adapterclasses ofsame actual loadedstock contributions; whole bodies unexpanded',K.OLIGO)
 other_material=g.product('tag_other_material','Other copresent source stockcomponents; loadingoccupancy/freecompanions unassembled',K.MOLECULAR)
 g.react('tag_components','Same selectedloaded-stock inventory, nofresh enzyme/adapter additions',S.GENERIC,((selected,R.TEMPLATE),(other,R.STRUCTURAL)),(enzyme,representative,adapters,other_material),changes_molecule=False)
 b.op('tag_array','Source distribute/allocate nuclei andloadedstocks BEFORE48wellreaction array (10,000nuclei/well,12barcodes×4wells);55C30min',ins=((enzyme,R.CATALYST),(representative,R.STRUCTURAL),(adapters,R.STRUCTURAL),(other_material,R.STRUCTURAL),b.cond('tag_conditions','Actual source Tn5 buffer;48well preparation, individual dose allocations unassembled')),product='Actual48reaction array collection after tagmentation; not an asserted single physical tube or one universal barcode')
 b.op('stop_tag','Source15mMEDTA terminates tagmentation',S.GENERIC,(b.stock('edta','Actual source EDTA'),))
 b.clean('pool_tag','Actual source nuclei pool/20umfilter/1000g10min4C recovery; array-to-pool handling compound')
 b.op('barcode1','Source distribute/allocate40ulperBCPlate01well BEFORE37C30minreaction array; actual member inventory unexpanded',ins=barcode_supplies(g,'ligation1','Source BCPlate01 barcodes/LinkerR01/T4ligase/BSA/buffer')+(b.cond('ligation1_conditions','Published37C30min; occupied nicked strands unassembled'),))
 block_stock=g.input('block1','Actualsource10ulR01blocking solution stock',K.MOLECULAR)
 blocker=g.product('block1_oligo','Copresent BlockerR01 component ofsameblocking solution',K.OLIGO,state=single(BLOCK1,'BlockerR01'))
 buffer=g.product('block1_buffer','CopresentsourceT4ligasebuffer/water component ofsameblocking solution',K.CONDITION)
 g.react('block1_components','ActualR01blocking-stock components, notfreshindependent additions',S.GENERIC,((block_stock,R.TEMPLATE),),(blocker,buffer),changes_molecule=False)
 b.op('block_round1','Source10ulR01blocking solution30min',S.HYBRIDIZATION,((blocker,R.STRUCTURAL),(buffer,R.CONDITION),b.cond('block1_conditions','Source30min')))
 b.clean('pool1','Actual completedBCPlate01reaction array→physical nuclei pool and1000g10min4C recovery before BCPlate02')
 b.op('barcode2','Source distribute/allocate toBCPlate02 BEFOREsecond37C30minreaction array; actual member inventory unexpanded',ins=barcode_supplies(g,'ligation2','Source BCPlate02 barcodes/LinkerR02/T4ligase/BSA/buffer')+(b.cond('ligation2_conditions','Source similar ligation schedule; full plate assemblies unexpanded'),))
 stopstock=g.input('termination','Actual source BlockerR02/EDTA termination solution',K.MOLECULAR)
 stopol=g.product('termination_oligo','Copresent BlockerR02 component ofsameactual solution',K.OLIGO,state=single(BLOCK2,'BlockerR02'))
 stopchem=g.product('termination_edta','Copresentactual EDTA/water component ofsame termination solution',K.MOLECULAR)
 g.react('termination_components','Actual termination-stock components, notextra freeblocker/EDTA additions',S.GENERIC,((stopstock,R.TEMPLATE),),(stopol,stopchem),changes_molecule=False)
 b.op('stop_barcode','Source BlockerR02/EDTA termination',S.GENERIC,((stopol,R.STRUCTURAL),(stopchem,R.STRUCTURAL)))
 b.clean('nuclei_count','Actual source PBS resuspension/count before sublibrary selection')
 b.select('sublibrary','Actual~2,000nuclei sublibrary in35ulPBS; other sublibraries outside this trace')
 b.op('lyse','Source55C2h sublibrary proteinaseK lysis',S.DEPROTEINIZATION,(b.stock('lysis_supplies','Actual source5ul4MNaCl+5ul10percentSDS supplies'),b.stock('lysis_protease','Actual separate5ul10mgmlProteinaseK',K.CATALYST,R.CATALYST),b.cond('lysis_conditions','Published55C2h850rpm')))
 b.clean('dna_clean','Source1xAMPure,30ulwater elution')
 b.op('denature1','Source72C5minthenice3min',S.THERMAL_TREATMENT,(b.cond('denature1_conditions','Source72C/ice sequence'),))
 b.op('exoi1','Source ExonucleaseI removes excessLinkerR02,37C60min',S.DIGESTION,(b.stock('exoi','Actual1ulExonucleaseI',K.CATALYST,R.CATALYST),b.cond('exoi_conditions','Actual1xExonucleaseIbuffer37C60min')))
 b.clean('exo_clean','Source1.8xAMPure')
 b.op('denature2','Source95C1minthenice',S.THERMAL_TREATMENT,(b.cond('denature2_conditions','Source95C1minthenice'),))
 b.op('bridge_ligate','Source3prime-end bridge ligation20C16h',S.LIGATION,(b.stock('bridge_adapter','Actual source BridgeAdapter10uM;5Ptop/N6bottom3NH2',K.ADAPTER,state=bridge()),b.stock('bridge_ligase','Actual30UT4DNA ligase',K.CATALYST,R.CATALYST),b.stock('bridge_peg','Actual7.5percentPEG8000'),b.cond('bridge_conditions','Source30ultotal in1xT4ligasebuffer20C16h')),geometry=AdapterLigationGeometry.known('splinted',bridge(),statement='Source-known local paired bridge/N6capture end, top5P/bottom3amine. Occupied genomic target, randomN6 pairing and sealed junction unassembled.'))
 b.clean('bridge_clean','Source1xAMPure,26ulwater')
 b.op('add_spike','Source10pgadaptor-ligatedT1M/T2H/T3MH control mixture',S.GENERIC,(b.stock('spike','Actual methylated/hydroxymethylated adaptor-ligatedDNA control mixture; source whole controls unexpanded',K.MOLECULAR,R.SPIKE_IN),))
 b.clean('preoxidation_beads','Source1xAMPure preoxidation purification')
 b.clean('p6','SourceMicroBioSpinP6 cleanup')
 b.op('denature_chemical','Source.05MNaOH37C30minthenice5min',S.CHEMICAL_CONVERSION,(b.stock('naoh','ActualsourceNaOH'),b.cond('denature_chemical_conditions','Source37C30minthenice5min')))
 b.op('oxidize_hmc','SourceK2RuO4 oxidation onice1h',S.CHEMICAL_CONVERSION,(b.stock('ruthenate','Actual K2RuO4 oxidant preparedfromperruthenate/NaOH25C2d; selected fresh1.5ul1x stock'),b.cond('oxidation_conditions','Source28.5ulDNA+1.5uloxidant,onice1h')),chemical={'source_relation':'5hmC oxidation to5fC class; no universal conversion/yield or whole sequence assigned'})
 b.clean('oxide_p6','SourceMicroBioSpinP6 after oxidation')
 b.op('label_hmc','Source150mMmalononitrile37C20h',S.CHEMICAL_CONVERSION,(b.stock('malononitrile','Actual150mMmalononitrile'),b.cond('malono_conditions','Source35ultotal10mMTrispH7/37C20h850rpm')),chemical={'source_relation':'Malononitrile labels5fC derivedfrom5hmC; later polymerase copying reports signal, not direct C-to-T DNA substitution atthis step'})
 b.op('record_hmc','SourceKAPA2GRobustHS indicator-primer extension,16cycles65C/68C',S.EXTENSION,(b.stock('indicator_primer','Actual two5caC indicator primer',K.OLIGO,state=indicator(),chemical_state={'modified_bases':({'zero_based_index':len(INDICATOR)-6,'modification':'5caC'},{'zero_based_index':len(INDICATOR)-2,'modification':'5caC'}),'source_ordering':RAW_ORDERINGS['P7-indicator primer']}),b.stock('extension_nt','Actual1ul10mMdNTPs'),b.stock('robust','Actual.8ulKAPA2GRobustHS',K.CATALYST,R.CATALYST),b.cond('extension_conditions','Source10ul5xKAPAbufferGC,95C5min/16×65C60s68C8min/72C10min')),chemical={'source_relation':'New complementary DNA copy records5hmC-derived signal; original genomicstrand distinct fromcopied strand. Occupied product/quantitative readout unassembled'})
 b.op('exoi2','Source ExonucleaseI37C60min',S.DIGESTION,(b.stock('exoi2_enzyme','Actual fresh1ulExonucleaseI',K.CATALYST,R.CATALYST),b.cond('exoi2_conditions','Source3ul10xbuffer37C60min')))
 b.clean('extension_clean','Source1.8xAMPure,10ulwater')
 b.op('oxidize_mc','SourcemTET1CD oxidation37C80min',S.CHEMICAL_CONVERSION,(b.stock('tet_enzyme','Actual8uMmTET1CD, purified enzyme stock; upstream expression not modelled',K.CATALYST,R.CATALYST),b.stock('tet_cosubstrates','ActualalphaKG/ascorbate/ATP/DTT/Fe2+ molecularsupplies'),b.cond('tet_conditions','Source50mMHEPES8/100mMNaCl/37C80min')),chemical={'source_relation':'Source5mC oxidation class before pic-borane; occupied intermediates and conversion yield unassembled'})
 b.op('proteolyse','Source.8UQiagenProtease50C1h',S.DEPROTEINIZATION,(b.stock('protease','ActualsourceQiagenProtease',K.CATALYST,R.CATALYST),b.cond('protease_conditions','Source50C1h')))
 b.clean('tet_clean','Source1.8xAMPure')
 b.op('reduce_mc','Sourcepic-borane70C4h',S.CHEMICAL_CONVERSION,(b.stock('borane','Actual12.5ul3.26Mpic-borane/DMSO+2.5ul3MsodiumacetatepH5.2'),b.cond('borane_conditions','Source10ulDNA/70C4h')),chemical={'source_relation':'Source original-strand5mC-derived products and primed5caC indicator treatment; no direct literal DNA C-to-T substitution or predicted quantitative completeness'})
 b.clean('borane_clean','SourceZymoOligoCleanConcentrator OR P30; actual choiceunassigned')
 # Known actual commercial master mix: its copresent functions, not fresh separate doses.
 supplies=g.input('pcr_supplies','Actual25ul2xKAPAU+HiFi mix plus2.5ul10uMUniversal andIndex primer contributions; original stockpackaging known for mix only',K.MOLECULAR)
 f=g.product('pcr_forward','SameactualsourceUniversal primer',K.OLIGO,state=single(UNIVERSAL,'Universal'))
 r=g.product('pcr_reverse','SameactualsourceP7indexprimer, selected6nt indexvalueunassigned',K.OLIGO,state=index_primer())
 enzyme=g.product('pcr_enzyme','SameKAPAU+HiFi master mix polymerase component',K.CATALYST);nt=g.product('pcr_nt','SameKAPAU+HiFi master mix nucleotide components',K.MOLECULAR);other=g.product('pcr_other','Otheractualmastermixsupplies',K.MOLECULAR)
 g.react('pcr_components','ActualPCRsupplyinventory, notextraindividualfreshdoses',S.GENERIC,((supplies,R.TEMPLATE),),(f,r,enzyme,nt,other),changes_molecule=False)
 b.op('index_pcr','Source4–6cycles finalindexPCR',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST),(nt,R.STRUCTURAL),(other,R.STRUCTURAL),b.cond('pcr_conditions','Source98C10s;4–6×98C10s68C15s72C1min;72C10min')),geometry=PCRGeometry.unavailable(Scene(),reason='SourcePCRprimerordersknown, completeconvertedgenomiclibrary/landinggeometry unassembled.'))
 b.clean('final_clean','Source1.8xAMPure,20ulwater')
 b.select('loading','Actual selected library contribution toMGISEQ2000 orNovaSeq; loadingvolume/pooling scheduleunassigned')
 library=g.product('library','Actual source selected PCR-DNA sequencing library; whole occupied constructs unassembled',K.LIBRARY,chemical_state={'source_processing_umi':feature('simple.umi','umi','random',note='Published processing defines first10basesRead2 asUMI; physicalinstallation/wholemember/primerlanding unassembled, notasserted specificallytoR02inthisselectedinventory')})
 g.react('readiness','Library readiness interface after actual sourcePCR/cleanup, nota new molecular transformation',S.GENERIC,((b.x,R.TEMPLATE),),(library,),changes_molecule=False)
 g.sequence(library,SequencingPlan(acquisition='unspecified',limitation='Published MGI200+6+100 orNovaSeq150+8+8+150. Exactrun/samplechoice, primerlandings, indexvalues andNovaSeqIndex2content unassigned; no computedreadmask.'))
 g.declare_reviewed_scope(BIO_SCOPE,AUTHORITY);g.declare_reviewed_scope(BARCODE_SCOPE,AUTHORITY);g.declare_reviewed_scope(FINISH_SCOPE,AUTHORITY)
 return g

def selftest():
 from checks import Check
 c=Check();g=protocol_graph()
 c('pMENTs modeled polymer includes additional terminal ddC nucleotide',''.join(x.top for x in g.materials['tn5_bottom'].state.strands['pMENTs'].segs),'CTGTCTCTTATACACATCTC')
 c('bridge whole fixed source top ordering preserved',BRIDGE_TOP,'CGTCAGCTGTCTCTTATACACATCTGACGCTGCCGACGA')
 c('indicator source modified ordering fixed bases preserved',INDICATOR,'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTGACG')
 c('tagging components derive from selected loaded stock',g.steps['tag_components'].inputs[0].material,'tn5_selected')
 c('tagging does not add newfree Tn5 enzyme',g.steps['tag_array'].inputs[1].material,'tag_enzyme')
 c('factory stock has separate untransferred remainder',g.terminals['tn5_remainder'].disposition,Disposition.UNRESOLVED)
 c('selected sublibrary leaves other nuclei outside trace',g.terminals['sublibrary_remainder'].disposition,Disposition.UNRESOLVED)
 c('bridge ligation after sourceExoIcleanup',g.steps['bridge_ligate'].inputs[0].material,'denature2_product')
 c('hmC recorder follows malononitrile',g.steps['record_hmc'].inputs[0].material,'label_hmc_product')
 c('5mC treatment follows indicator copy/ExoIcleanup',g.steps['oxidize_mc'].inputs[0].material,'extension_clean_kept')
 c('known linker carried fromactualbarcode supplies',g.steps['barcode1'].inputs[2].material,'ligation1_linker')
 c('finalPCR uses produced sourceprimer inventory',g.steps['index_pcr'].inputs[1].material,'pcr_forward')
 c('loading leaves untransferred finalDNA stock',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
 c('wholePCRgeometry remains explicitunavailable',g.steps['index_pcr'].geometry.availability,'unavailable')
 from checks import Source,have
 from storage import sources_root
 source=sources_root()/'simple-seq__10.1038+s41587-024-02148-9'/'s41587-024-02148-9_41587_2024_2148_MOESM2_ESM.xlsx.txt'
 if have(c,Source(source,'Published SupplementaryTable1 original: https://doi.org/10.1038/s41587-024-02148-9; tools/doctext.py produces faithful selected-table twin'),label='Produced oligo whole-source transcription'):
  # Only selected Table1_1, whose raw original ZIP/XML was independently verified.
  # Do not scan the unselected 29MB result/barcode sheets or restate model constants.
  rows={}
  with source.open() as handle:
   for line in handle:
    if line.startswith('## sheet ') and 'Supplementary_Table1_1' not in line:break
    fields=line.rstrip('\n').split('\t')
    if len(fields)>1:rows[fields[0]]=fields[1]
  for target,name in (('ligation1_linker','Linker-R01'),('ligation2_linker','Linker-R02'),('pcr_forward','N5 Universal primer')):
   state=g.materials[target].state
   modeled=''.join(segment.top for strand in state.strands.values() for segment in strand.segs)
   c('Produced whole oligo matches original '+name,modeled,rows.get(name,''))
 c.report()
if __name__=='__main__':selftest()

PAGE_FILE='simple_seq.html'
def render_page():
 from page import head,info,caveat,table
 from chemdraw import panel
 from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
 from chemistry_validation import check_graph,Status
 g=protocol_graph();g.require_publishable()
 start='<div class="wrap" data-protocol-complete="true"><h1>SIMPLE-seq</h1><p>Measure 5mC and 5hmC from the same nucleus using combinatorial barcoding and sequential chemical recording.</p>'
 start+=info('Source: <a href="https://doi.org/10.1038/s41587-024-02148-9">Bai et al., Nature Biotechnology</a>, published online9February2024. Published Methods and original SupplementaryTable1 govern this page. The authors’ bench-protocol and GEO alternatives are not substituted.')
 start+=caveat('This is a source-qualified physical preparation with selected known molecular interfaces. Individual barcode-pool members, occupied genomic molecules, conversion fractions, complete ligation/PCR products and sequencing-primer landings remain unassembled. The48-well tagmentation and subsequent barcode plates are explicit compound arrays, not one physical reaction or a fabricated single-barcode experiment.')
 start+='<h2>Preparation</h2>'+table(('Stage','Published route'),(
 ('Barcode nuclei','Fix/isolate nuclei; load12barcoded Tn5stocks; tagment48wells (fourwells/barcode), pool, then perform two ligation-barcode rounds with blocking. Select~2,000nuclei per sublibrary.'),
 ('Tag genomic DNA','ProteinaseK lysis and cleanup; ExonucleaseI workup; denature and ligate the source5′-phosphorylated bridge with a randomN6 capture overhang and3′-amine-blocked partner.'),
 ('Record5hmC','Add controlDNA; NaOH denature; K2RuO4 oxidation and malononitrile labelling. KAPA2GRobustHS extends the two-5caC indicator primer, producing complementary copied DNA.'),
 ('Record5mC','Treat with mTET1CD, protease and pic-borane. Original-strand and copied-strand signals remain distinct; chemistry is not implemented as universal literal C-to-T replacement.'),
 ('Finish','KAPAU+HiFi indexPCR4–6cycles,1.8×AMPure cleanup and source20µl elution. A selected library contribution goes to sequencing; untransferred stock remains explicit.'),))
 start+='<h2>Known local oligo interfaces</h2>'+panel(transposon().rows(),caption='Source representative DNA01/pMENTs pair:19bpME stem; the additional terminalddC is retained and explicitly outside asserted pairing.')+panel(bridge().rows(),caption='Source bridge duplex: paired fixedcore, top5′P and bottomN6/3′NH2. Genomic capture occupancy and sealed genomic junction are unassembled.')
 start+=table(('Component','Source ordering'),tuple(RAW_ORDERINGS.items()))
 start+='<h2>Read setup and source differences</h2>'+table(('Source','Declaration'),(
 ('Published MGI','MGISEQ2000: Read1 200,Index1 6,Read2 100.'),('Published Illumina','NovaSeq: Read1 150,Index1 8,Index2 8,Read2 150.'),('Source processing','First10basesRead2 are UMIs; three barcode components jointly identify cells. Whole member installation/occupied layout remains unassembled.'),('Unresolved assignment','Exact run/sample choice, index values, primer landings and NovaSeqIndex2 contents are unassigned. No computed complete read mask is claimed.'),))
 start+=caveat('Published Methods print5Phos on barcoded Tn5 DNA, while the workbook barcode rows omit that modification. The model uses the explicit Methods5′P and retains the difference. Methods specify KAPA2GRobustHS,48wells and~2,000nuclei/sub-library; the authors’ separately dated bench protocol uses KAPAHiFi,96wells and5,000–10,000nuclei. These documents are not treated as equivalent versions. Historical GEO protocols give MGI150+8+150; the published layouts above govern this page, but individual deposited-run assignment remains unresolved. The source6ntP7index versus8index cycles and index-free Universal primer versus NovaSeqIndex2 are not silently repaired.')
 claims=check_graph(g)
 start+='<h2>Validation scope</h2>'+table(('Molecular claim','Count'),tuple((x.value,str(sum(c.status==x for c in claims)))for x in Status))
 return assemble_page(g,[head('SIMPLE-seq — joint5mC/5hmC preparation'),start,'<details><summary>Complete preparation, inventories and unresolved transitions</summary>',WorkflowSection('Full source preparation and physical fates'),PlanLimitationSection('library'),'</details><p><a href="01_simple-seq.html">Research note and source differences</a></p></div>'])
