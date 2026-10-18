"""Published plate-based Prox-seq: actual ssDNA reporters, no assumed RT."""
from chemdraw import Construct,Scene,Segment,feature,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import seqprimers as sp
TITLE='Prox-seq — published plate-based reporter preparation'
NOTES='01_prox-seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-022-01684-z">Vistain etal.,NatureMethods2022</a>; publishedMOESM1 plateMethods and MOESM3Tables7–10.'
SUMMARY='A/B antibody-linked DNA probes ligate on a connector. In the selected non-PFA plate experiment, single-cell lysate directly feeds two DNA PCRs; the source did not measure mRNA in that plate experiment.'
CAVEAT='This graph selects the source-positive non-PFA plate implementation rather than mixing Drop-seq,10x and plate chemistry. Other published implementations and their pairedRNA readouts remain distinct; their native RT/cDNA/physicalfractionation routes are deferred here. SourceTable10 HHHHACHHHHACHHHN is a constraineddegenerate probeUMI (H=A/C/T), with the printed fixedAC interruptions. Probe azide/antibody linkages remain on originals, not their new PCR copies. Full nativefirst-cycle/endgeometry remains bounded; no final expectedduplex is used as template. Optionallysate reservation does not establish an observed expressionlibrary.'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_ENDING='Selectedplate: NextSeqmid-outputv2.5,40%PhiX, customSmartPLARead1/i5/i7 primers;75/8/8, noRead2. Drop-seq20/8/85 and10xPLA28/8/75 are separate published implementations.'
READ_SPECS=(sp.read_spec('Read 1',75,'probeAUMI and pairedproteinbarcodes; customSmartPLA_Read1','extract'),sp.read_spec('Read 2','not acquired','Sourceplate acquisition is single-ended','undetermined'),sp.read_spec('Index 1 (i7)',8,'singlecell librarydualindex','demultiplex'),sp.read_spec('Index 2 (i5)',8,'singlecell librarydualindex','demultiplex'))
PROBE_A_5_SITE = 'CGCATTGCATCGTCTCGTGGGCTCGGC'
PROBE_A_UMI = 'HHHHACHHHHACHHHN'
PROBE_A_SPACER = 'GCAG'
PROBE_A_LINKER = 'GATCG'
PROBE_A_ARM = 'CTAAATCGTG'
PROBE_B_ARM = 'TCGTGTCGTG'
PROBE_B_LINKER = 'TCTAAAGTCC'
PROBE_B_SPACER = 'AC'
PROBE_B_POLYA = 'A' * 29
BARCODE_LEN = 8
CONNECTOR = 'TTTCACGACACGACACGATTTAGGTC'
PROBE_A_LEN = len(PROBE_A_5_SITE) + len(PROBE_A_UMI) + len(PROBE_A_SPACER) + BARCODE_LEN + len(PROBE_A_LINKER) + len(PROBE_A_ARM)
PROBE_B_LEN = len(PROBE_B_ARM) + len(PROBE_B_LINKER) + BARCODE_LEN + len(PROBE_B_SPACER) + len(PROBE_B_POLYA)
PRODUCT_LEN = PROBE_A_LEN + PROBE_B_LEN
PROBE_A_TABLE = '/5AzideN/' + PROBE_A_5_SITE + PROBE_A_UMI + PROBE_A_SPACER + '[8-nt barcode]' + PROBE_A_LINKER + PROBE_A_ARM
PROBE_B_TABLE = '/5Phos/' + PROBE_B_ARM + PROBE_B_LINKER + '[8-nt barcode]' + PROBE_B_SPACER + PROBE_B_POLYA + '/3AzideN/'
BARCODE_A = feature('protein_barcode_a', 'feature_barcode', 'whitelist',whitelist='publishedTable10antibodybarcodefamily')
BARCODE_B = feature('protein_barcode_b', 'feature_barcode', 'whitelist',whitelist='publishedTable10antibodybarcodefamily')
PLA_UMI = feature('pla_umi', 'umi', 'degenerate', note='SourceHHHHACHHHHACHHHN constrainedwindow: H=A/C/T, N=A/C/G/T; fixedAC interruptions retained')
def probe_a_parts():
    return (Segment('probe-A U_fwd site', PROBE_A_5_SITE), Segment('PLA UMI', PROBE_A_UMI, 'umi', placeholder=True, feature=PLA_UMI), Segment('probe-A spacer', PROBE_A_SPACER), Segment('antibody A barcode', 'N' * BARCODE_LEN, 'cbc', placeholder=True, feature=BARCODE_A), Segment('probe-A linker', PROBE_A_LINKER), Segment('A connector arm', PROBE_A_ARM))
def probe_b_parts():
    return (Segment('B connector arm', PROBE_B_ARM), Segment('probe-B linker', PROBE_B_LINKER), Segment('antibody B barcode', 'N' * BARCODE_LEN, 'cbc', placeholder=True, feature=BARCODE_B), Segment('probe-B spacer', PROBE_B_SPACER), Segment('poly(A)', PROBE_B_POLYA))
def connector_geometry():
    """Table 10 probe ends on the published connector; Scene checks all 20 pairs."""
    sc = Scene()
    sc.strand('probe A end', [Segment('probe-A linker', PROBE_A_LINKER), Segment('A connector arm', PROBE_A_ARM)], label='probe A 3′ end (3′-OH)')
    sc.anneal('connector', [Segment('connector 5′ overhang', CONNECTOR[:3]), Segment('connector B site', CONNECTOR[3:13]), Segment('connector A site', CONNECTOR[13:23]), Segment('connector 3′ overhang', CONNECTOR[23:])], to='probe A end', pair=('connector A site', 'A connector arm'), label='universal connector (26 nt)', unpaired=('connector 5′ overhang', 'connector 3′ overhang'))
    sc.anneal('probe B end', [Segment('B connector arm', PROBE_B_ARM), Segment('probe-B linker', PROBE_B_LINKER)], to='connector', pair=('B connector arm', 'connector B site'), label='probe B 5′ end', mod5='p', above=True, unpaired=('probe-B linker',))
    sc.ligation_end('probe B end', 'B connector arm', kind='splinted', polarity='5-prime', owner='adapter', text='5′-phosphate nick partner')
    return AdapterLigationGeometry.known('splinted', sc, statement='Supplementary Table 10 probe A ends in CTAAATCGTG-3′ and probe B begins 5′-phosphate-TCGTGTCGTG; the 26-nt connector pairs with both 10-nt arms (20 bases, as the article states) and leaves 3-nt single-stranded ends. The arm assignment is computed from the printed sequences; the paper does not draw it.')

RT_HANDLE='AAGCAGTGGTATCAACGCAGAGTAC'
U_FWD='GCATCGTCTCGTGGGCTC'
P5_ARM='TATTTGGCATCGTCTCGTGGGCTCGG'
P7_ARM='GCCGAAGCAGTGGTATCAACGCAGAGT'
I5=feature('prox.platei5','sample_index','whitelist',whitelist='publishedTable9celllibraryindices')
I7=feature('prox.platei7','sample_index','whitelist',whitelist='publishedTable9celllibraryindices')
# SourceTable8i5family/Table9N701; no claim everywell uses this pair.
P5_PREFIX='AATGATACGGCGACCACCGAGATCTACAC'
P7_PREFIX='CAAGCAGAAGACGGCATACGAGAT'
P5_INDEX='[i5]'
P7_INDEX='TCGCCTTA'
def index_parts(p5):return [Segment('sourceP5prefix' if p5 else 'sourceP7prefix',P5_PREFIX if p5 else P7_PREFIX),role('sourcei5member/physicalwidthunassigned',molecular_feature=I5) if p5 else Segment('sourcecelllibraryindex',P7_INDEX,feature=I7),Segment('sourceannealingarm',P5_ARM if p5 else P7_ARM)]
def oligo_state(name,parts,mod5='',mod3=''):
 sc=Scene();sc.strand(name,parts,mod5=mod5,mod3=mod3,label='sourceoligo');return sc

def mix(stage):
 sc=Scene()
 if stage in ('cell','bound','hybridized','joined','lysate'):
  sc.strand('RNA',[role('originalcellularRNA; physicalU')],label='independentRNA population')
  sc.strand('genome',[role('originalgenomicDNA')],label='independentgenome')
  sc.strand('proteins',[role('cellular/antibodyproteins; attachmentgeometryunassigned')],label='independentprotein components; no DNA–protein concatenation')
 if stage in ('conjugate','bound','hybridized'):
  sc.strand('A',probe_a_parts(),mod5='antibody conjugate',label='originalA:5′azide clicked toDBCOantibody')
  sc.strand('B',probe_b_parts(),mod5='p',mod3='antibody conjugate',label='originalB:3′azide clicked toDBCOantibody')
 if stage in ('joined','lysate'):
  sc.strand('reporter',[*probe_a_parts(),*probe_b_parts()],label='originalcovalentA→B singleDNA reporter; antibodycouplingretained')
  sc.junction('reporter','A connector arm','B connector arm',text='sourceconnector-splinted nicksealed; connectornotinserted')
 if stage in ('preamp','indexed','pooled'):
  parts=[role('nativePCRfootprint; fullassemblybounded'),*probe_a_parts()[1:],*probe_b_parts(),role('copiedoligodTGT handle')]
  if stage in ('indexed','pooled'):parts += [role('installedsampleindexi5; physicalwidthunassigned',molecular_feature=I5),role('installedsampleindexi7',length_bp=8,molecular_feature=I7)]
  sc.strand('daughter',parts,label='new normalDNA reporter copies; originalazide/antibody notcopied')
 if stage=='conjugate':sc.strand('proteins',[role('independentDBCOantibody components; DNAendcouplingnotrenderedascovalentproteinsequence')],label='antibodyconjugates')
 if stage not in ('cell','conjugate'):sc.strand('residuals',[role('possible originalreporters/unligatedprobes/genomicDNA/RNA/antibodies/oligos andearliercopies; laterphysicalrecovery unassigned')],label='independent residualorigins; nonamplification not physicaldiscard')
 return sc

def oligos():
 return [f'<b>ProbeA</b>:5′-{PROBE_A_TABLE}-3′',f'<b>ProbeB</b>:5′-{PROBE_B_TABLE}-3′',f'<b>Connector</b>:5′-{CONNECTOR}-3′',f'<b>OligodTGT</b>:5′-{RT_HANDLE+"T"*30+"GT"}-3′',f'<b>U_fwd</b>:5′-{U_FWD}-3′',f'<b>SmartPLA_P5_S5XX</b>:5′-{P5_PREFIX+P5_INDEX+P5_ARM}-3′',f'<b>SmartPLA_P7_N701</b>:5′-{P7_PREFIX+P7_INDEX+P7_ARM}-3′']

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def clean(i,l,x):
  sel=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+' retainedfraction',M.MOLECULAR,state=x.state);w=g.product(i+'_waste','actual nonretainedsupernatant/washes',M.WASTE);g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));g.discard(w,'Actual removedfraction; complete molecularclass depletion not presumed.');return y
 cell=g.input('cell','non-PFAfixed Jurkat/Raji sourceplate cells',M.BIOLOGICAL,state=mix('cell'),primary=True)
 antibodies=inp('antibodies','DBCO-conjugatedantibodies; sourceDBCO NHSreaction/50kDa cleanup; exactstoichiometryunassigned',M.MOLECULAR)
 conjugation_conditions=inp('conjugation_conditions','source4°Covernight')
 a=inp('probe_a','Table10probeA5′azideDNA',M.OLIGO,oligo_state('A',probe_a_parts(),mod5='azide'))
 b=inp('probe_b','Table10probeB5′P/3′azideDNA',M.OLIGO,oligo_state('B',probe_b_parts(),mod5='p',mod3='azide'))
 conjugates=g.product('conjugates','A/B DNAconjugatedantibodymixture',M.MOLECULAR,state=mix('conjugate'))
 g.react('conjugate','CoupleazideDNA toDBCOantibodies;4°Covernight',K.CHEMICAL_CONVERSION,((a,I.TEMPLATE),(b,I.TEMPLATE),(antibodies,I.STRUCTURAL),(conjugation_conditions,I.CONDITION)),(conjugates,))
 binding=inp('binding_conditions','2.5nMeachA/B;37°C60min',M.CONDITION)
 bound=g.product('bound','probe-boundcells',M.MOLECULAR,state=mix('bound'))
 g.react('bind','BindDNA-conjugatedantibodies',K.BINDING,((cell,I.TEMPLATE),(conjugates,I.STRUCTURAL),(binding,I.CONDITION)),(bound,))
 washed=clean('binding_washes','Threepellet/resuspend1%BSA/PBS washes',bound)
 connector=inp('connector','source26ntuniversalconnector',M.OLIGO,oligo_state('connector',[Segment('sourceconnector',CONNECTOR)]))
 hybrid=g.product('hybrid','connector-associatedproximalprobes',M.MOLECULAR,state=mix('hybridized'))
 g.react('hybridize','Hybridizetheconnector toA/B10ntarms',K.HYBRIDIZATION,((washed,I.TEMPLATE),(connector,I.STRUCTURAL)),(hybrid,))
 ligase=inp('ligase','130U/mlT4DNAligase,HEPES/Mg/rATP;37°C3h',M.CATALYST)
 joined=g.product('joined','covalentreporters oncells withresidualunligatedprobes',M.MOLECULAR,state=mix('joined'))
 g.react('ligate','SealA3′OH toB5′P onconnector',K.LIGATION,((hybrid,I.TEMPLATE),(connector,I.STRUCTURAL),(ligase,I.CATALYST)),(joined,),geometry=connector_geometry())
 joined=clean('post_ligation','Pelletandresuspend1%BSA/PBS+PI beforelivesort',joined)
 sort=inp('sort','SinglePI-negativecell per96well;0.1%Tritonlysis/oligodTprimers/dNTPs;72°C3min',M.CONDITION)
 lysate=g.product('lysate','actualsinglecelllysate; independentRNA/DNA/reporter populations',M.MOLECULAR,state=mix('lysate'))
 unallocated=g.product('unallocated','unsorted/nonselectedsourcecells',M.MOLECULAR)
 g.react('sort_lyse','Sortsinglecells andheatlysate',K.PARTITION,((joined,I.TEMPLATE),(sort,I.CONDITION)),(lysate,unallocated));g.mark_terminal(Terminal(unallocated.id,Disposition.UNRESOLVED,'Notselectedsourceallocation; not assumed physicallydestroyed.'))
 used=g.product('preamp_aliquot','2–4µl actualmixedlysate forPLApreamp',M.MOLECULAR,state=lysate.state)
 remainder=g.product('remaining_lysate','remainingmixedlysate; optionalSmartseq2notperformed inselectedexperiment',M.MOLECULAR,state=lysate.state)
 aliquot_rule=inp('aliquot_rule','source2–4µl mixedlysate')
 g.react('allocate','Takeactuallysate aliquot; no modalitypurification',K.SAMPLING,((lysate,I.TEMPLATE),(aliquot_rule,I.CONDITION)),(used,remainder));g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Paper saysplateexperiment didnotmeasuremRNA; possible optionalroute is not an observedlibrary.'))
 f=inp('u_fwd','0.1µMU_fwd',M.OLIGO,oligo_state('U_fwd',[Segment('U_fwd',U_FWD)]))
 r=inp('dtgt','0.1µMoligodTGT',M.OLIGO,oligo_state('oligodTGT',[Segment('handle',RT_HANDLE),Segment('dT30','T'*30),Segment('GT','GT')]))
 pol=inp('kapa','KAPAHiFiHotStartReadymix',M.CATALYST)
 pre=g.product('preamp','newreporterDNA pluspossibleoriginals',M.MOLECULAR,state=mix('preamp'))
 g.react('preamp_pcr','DirectDNApreamp:5cycles55°C+17cycles67°C;noRT inselectedplatepath',K.PCR,((used,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(pre,),geometry=PCRGeometry.unavailable(used.state,reason='Actualmixedlysatecontains single-strandcovalentA–Breporter; sequentialfirstcopy/fullprimergeometry bounded.'))
 pre=clean('preamp_cleanup','1.8×AMPurepreampcleanup',pre)
 i5=inp('index_i5','SmartPLA_P5_S5XX;Table8nativei5family;member/physicalwidthunassigned',M.OLIGO,oligo_state('i5',index_parts(True)))
 i7=inp('index_i7','SmartPLA_P7_N7XX;Table9N701example',M.OLIGO,oligo_state('i7',index_parts(False)))
 indexed=g.product('indexed','newduallyindexedPLAcopies pluspossibleoriginals',M.MOLECULAR,state=mix('indexed'))
 g.react('index_pcr','2µlpreampinput;12cycles67°C indexedPCR',K.PCR,((pre,I.TEMPLATE),(i5,I.STRUCTURAL),(i7,I.STRUCTURAL),(pol,I.CATALYST)),(indexed,),geometry=PCRGeometry.unavailable(pre.state,reason='Actualrecoveredpreamp molecules; no finalendpointtemplate shortcut.'))
 indexed=clean('indexed_cleanup','1.8×AMPureindexedcleanup',indexed)
 pool=inp('other_indexed_wells','independentotherwellsnormalizedequimolar',M.MOLECULAR,mix('indexed'))
 library=g.product('library','physicalequimolarpooledplatePLAlibrary',M.LIBRARY,state=mix('pooled'))
 g.react('pool','Normalizeandpoolequimolarwelllibraries',K.POOL,((indexed,I.TEMPLATE),(pool,I.TEMPLATE)),(library,),changes_molecule=False)
 g.sequence(library,SequencingPlan(read_lengths={'Read 1':75,'Index 1 (i7)':8,'Index 2 (i5)':8},acquired_roles=('Read 1','Index 1 (i7)','Index 2 (i5)'),limitation=SEQUENCING_ENDING+' Fullnativeendassembly/customprimerlandings bounded.'))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('aliquot_rule','conjugation_conditions','cell','antibodies','probe_a','probe_b','conjugate','binding_conditions','bind','binding_washes_selection','binding_washes','connector','hybridize','ligase','ligate','post_ligation_selection','post_ligation','sort','sort_lyse','allocate','u_fwd','dtgt','kapa','preamp_pcr','preamp_cleanup_selection','preamp_cleanup','index_i5','index_i7','index_pcr','indexed_cleanup_selection','indexed_cleanup','other_indexed_wells','pool','library')
