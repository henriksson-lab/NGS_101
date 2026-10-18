"""icSHAPE2016 native markedRNA, interruptedssDNA primer and circularPCR input."""
from chemdraw import Construct,Scene,Segment,MolecularState,feature,circle_rows
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
import seqprimers as sp
TITLE='icSHAPE — source-native clickSHAPE library'
NOTES='01_icshape.html'
SOURCE='<a href="https://doi.org/10.1038/nprot.2016.011">Flynn et al.,NatureProtocols2016</a>; published author-manuscript PMC4896316.'
SUMMARY='Selected in-vivo NAI-N3 branch: RNA2′O adducts are clicked tobiotin; RT occurs before streptavidin enrichment. A blockedadapter and interrupted5′phosphorylated RTprimer yield actual circularssDNA PCR templates.'
CAVEAT='NAI-N3-treated branch selected; DMSO controls use a distinctbiotinylatedadapter, not the ddC adapter here. DDDNN–AACC–NNNN is a physicallydiscontiguous9-base UMI around a4-base samplebarcode; D meansA/G/T. TwoiSp18 moieties are nonnucleotide spacers, not18DNA bases. Native circlefold/PCRfirstcycle/fullfinalassembly remain unassigned. SourceStep49 says1µl1µM P3 in50µl yet20mM final; dilution implies20nM. No correctedsourceconcentration is assigned.'
ADAPTER='AGATCGGAAGAGCGGTTCAG'
RTHANDLE='AGATCGGAAGAGCGTCGTGGA'
RTTAIL='TACTGAACCGC'
RT_ORDER='/5phos/DDDNNAACCNNNN'+RTHANDLE+'/iSp18/GGATCC/iSp18/'+RTTAIL
P5='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
P3='CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
SAMPLE=feature('icshape.sample','sample_index','fixed',note='AACC experimentalbarcode physicallybetween two randomUMIcomponents')
UMI5=feature('icshape.umi5','umi','degenerate',group='icshape.umi',part='1',note='DDDNN first5positions ofsame9ntUMI; D=A/G/T, notC')
UMI4=feature('icshape.umi4','umi','random',group='icshape.umi',part='2',note='NNNN final4positions ofsame9ntUMI; concatenatewithDDDNN, not a secondindependentUMI')
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
FINAL_LIBRARY=Construct([role('native P5PCRend'),role('13ntinlinebarcodeUMIregion; physicalDDDNN/AACC/NNNNprimerarchitecture'),role('RT-stop-derivedDNAinsert'),role('native P3PCRend')],name='icSHAPE native libraryreference')
FINAL_CAPTION='Native componentreference; actual circularssDNA, not this endpoint, feedsPCR.'
SEQUENCING_ENDING='SourceStep68:75-base single-end sufficient; omitstandard6/8ntindexreads because IDsareinline. FirstmappedRTstop is computationallyshifted toinferadduct; that coordinate is not an insertedbase.'
READ_SPECS=(sp.read_spec('Read 1',75,'inline13nt contains4experimentalbarcode+9degenerate positions; thenRT-stopDNA; fullnative runprimergeometry unassigned','extract_then_include'),sp.ReadContent('Index 1 (i7)','not acquired','SourceStep68 instructs omission ofstandardindexread','No — not acquired'),sp.ReadContent('Index 2 (i5)','not acquired','SourceStep68 instructs omission ofstandardindexread','No — not acquired'))
def rt_parts():
 return [Segment('UMIfirst5','DDDNN',placeholder=True,feature=UMI5,note='physicalDNA ambiguityclasses; D=A/G/T'),Segment('sampleAACC','AACC',feature=SAMPLE),Segment('UMIfinal4','NNNN',feature=UMI4),Segment('RThandle',RTHANDLE),role('firstiSp18 nonnucleotide spacer'),Segment('GGATCC','GGATCC'),role('secondiSp18 nonnucleotide spacer'),Segment('annealingRTtail',RTTAIL)]
def rt_primer():
 sc=Scene();sc.strand('RTprimer',rt_parts(),mod5='p',label='actual interruptedDNA primer; two internalchemical spacers');return sc

def adapter_state():
 sc=Scene();sc.strand('adapter',[Segment('sourceadapterDNA',ADAPTER),Segment('physicalterminalddC','C',note='source/3ddC/ addsphysicaldeoxyC; no3′OH')],mod5='App',mod3='ddC-blocked',label='NAI-N3branchadapter');sc.ligation_end('adapter','sourceadapterDNA',kind='single-strand',polarity='5-prime',owner='adapter',text='activated5′App acceptsRNA3′OH');return sc

def native_state(stage):
 sc=Scene()
 if stage in ('raw','modified','clicked','fragmented','repaired','ligated'):
  r=role('original cellularRNA; physicalU; '+('NAI-N3 2′Oadducts remain nucleotide substituents' if stage!='raw' else 'nativefold/modifications unassigned'))
  parts=[r]
  if stage=='ligated':parts += [Segment('ligatedadapter',ADAPTER),Segment('terminalddC','C',note='originaladapter ddC physicalC;3′blocked')]
  sc.strand('RNA',parts,mod3='ddC' if stage=='ligated' else '',label='originalRNA '+stage+'; biotin-clickedadducts' if stage in ('clicked','fragmented','repaired','ligated') else 'originalRNA '+stage)
  sc.strand('otherRNA',[role('independent unmodified/unclickedRNA and possible otherinputpolymers; enrichment not pureclass')])
 elif stage in ('rt','released'):
  sc.strand('cDNA',[*rt_parts(),role('DNAcopyofRNA ending atRTstop; exactstopcoordinateunknown')],mod5='p',label='actual native interruptedssDNA; sourceRNAadduct/biotin not copied')
  sc.strand('RNAlineage',[role('originalbiotinylatedRNA plus unmodifiedRNA' if stage=='rt' else 'possibleRNasehydrolysates/retainedRNA fragments; chemicaldigestion not physicalclearance')])
 elif stage=='circle':
  parts=[*rt_parts(),role('RT-stopDNAcopy')]
  rows=list(circle_rows(Construct(parts,name='native interruptedssDNAcircle'),closure_label='sourceCircLigase linksDNA3′OH toprimer5′P; originaliSp18spacers remain internal'))
  sc.strand('residual',[role('possible unreactedlinearssDNA, RNAhydrolysates, nativeoligos; no separationatclosure')]);rows.extend(sc.rows());return MolecularState('actual interruptedssDNAcircle plusresidualclasses',tuple(rows))
 elif stage=='library':
  sc.strand('daughter',[Segment('sourceP5primer',P5),role('nativePCR-selectedRTstop/inlineidentifierstrandorder; firstcycletopologyunassigned'),role('retainedDDDNNcompoundUMIcomponent',length_bp=5,molecular_feature=UMI5),role('retainedAACCsampleidentifier; orientationunassigned',length_bp=4,molecular_feature=SAMPLE),role('retainedNNNNcompoundUMIcomponent',length_bp=4,molecular_feature=UMI4),Segment('sourceP3primerderivedend','[P3-derivedcopiedend]',placeholder=True)],label='DNA PCRdaughters; originalRNAadduct/biotin not copied')
  sc.strand('nativeoriginals',[role('possible originalcircular/linearDNA, blockedadapter/RTprimer/spacer/offtargetclasses; recoveryunassigned')])
 sc.strand('residuals',[role('possible co-carried originalpolymers/oligofragments; classspecificrecovery unresolved')])
 return sc

def oligos():return [f'NAI-N3branch3′adapter—5′-/5rApp/{ADAPTER}/3ddC/-3′',f'RTprimer1—5′-{RT_ORDER}-3′',f'P5-Solexa—5′-{P5}-3′',f'P3-Solexa—5′-{P3}-3′']

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def prod(i,l,stage,kind=M.MOLECULAR):return g.product(i,l,kind,state=native_state(stage))
 def op(i,l,k,x,extra,state):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state);g.react(i,l,k,((x,I.TEMPLATE),*extra),(y,));return y
 def clean(i,l,x):
  s=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+'retainedmixture',x.kind,state=x.state);w=g.product(i+'_waste','physically nonretainedfraction',M.WASTE)
  g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(s,I.SELECTION)),(y,w));g.discard(w,'Actualsourcephysical fraction; classpurity/recovery not presumed.');return y
 raw=g.input('cellular_rna','RNA in independentlytreated cellaliquot; selectedNAI-N3branch',M.BIOLOGICAL,state=native_state('raw'),primary=True)
 nai=inp('nai_n3','100mMNAI-N3;37°C5min');x=op('acylate','Acylate flexibleRNA2′OH',K.CHEMICAL_CONVERSION,raw,((nai,I.CONDITION),),native_state('modified'))
 isolation=inp('rna_isolation','sourceRNAisolation; detailedSteps4–11notreconstructed')
 isolated=prod('isolated','isolated NAI-N3RNA mixture','modified')
 g.react('isolate_rna','NativeRNA isolation boundary beforepolyAselection',K.MAGIC,((x,I.TEMPLATE),(isolation,I.CONDITION)),(isolated,),reason='Knownsourceisolation precedesProcedure13. Detailedextraction chemistry4–11 remains outside this capped inspectedsection scope; no genericpurificationreagents assigned.')
 x=isolated
 for n in (1,2):x=clean('polya'+str(n),'PolyAPuristMAG round'+str(n)+'; sourcepolyAselection/wash',x)
 x=clean('polya_rneasy','Lyophilize/RNeasymini RNAcleanup',x)
 click=inp('dibo_biotin','1.85mMDIBO-biotin;37°C2h',M.MOLECULAR);cond=inp('click_condition','copper-free click; sourceStep19')
 x=op('click','Clickbiotin onRNA NAI-N3azides',K.CHEMICAL_CONVERSION,x,((click,I.STRUCTURAL),(cond,I.CONDITION)),native_state('clicked'));x=clean('click_cleanup','ZymoRNA C&C5;freeze-dry',x)
 frag=inp('fragmentation_condition','95°Cpreheat60s;fragmentationreagent30–60s');x=op('fragment','FragmentoriginalclickedRNA',K.FRAGMENTATION,x,((frag,I.CONDITION),),native_state('fragmented'))
 stop=inp('fragment_stop','sourcefragmentationstopsolution;ice');x=op('stop_fragment','Stop fragmentation',K.QUENCHING,x,((stop,I.CONDITION),),x.state);x=clean('fragment_cleanup','ZymoRNA C&C5;freeze-dry',x)
 repair=inp('end_repair','FastAP+T4PNK;37°C1h',M.CATALYST);x=op('repair','RepairRNA3′ends',K.END_REPAIR,x,((repair,I.CATALYST),),native_state('repaired'))
 ad=inp('adapter','sourceNAI-N3adapterApp plusphysical3′ddC',M.ADAPTER,adapter_state());lig=inp('t4_rnl1','sourceT4RNA ligase1 highconcentration;25°C3h',M.CATALYST)
 ligated=prod('ligated','RNA3′adapterligated originalRNA','ligated')
 g.react('ligate','Ligate RNA3′OH toadapter5′App',K.LIGATION,((x,I.TEMPLATE),(ad,I.STRUCTURAL),(lig,I.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known('single-strand',adapter_state(),statement='Procedure26 explicitlyT4RNA ligase1; activatedadapter5′App andphysical3′ddC blocking.'))
 x=clean('ligation_cleanup','ZymoRNA C&C5;freeze-dry',ligated);x=clean('rna_page','DenaturingPAGE isolateadapterligatedRNA >~40nt; SpinX/Amiconelution',x);x=clean('rna_page_column','ZymoRNA C&C5;freeze-dry',x)
 rp=inp('rt_primer','source5′P RTprimer1; DDDNN/AACC/NNNN and twoiSp18',M.OLIGO,rt_primer());rt=inp('ssiii','SuperScriptIII;25°C3min42°C5min52°C30min',M.CATALYST)
 x=op('rt','Reverse-transcribe RNA beforestreptavidinenrichment',K.REVERSE_TRANSCRIPTION,x,((rp,I.STRUCTURAL),(rt,I.CATALYST)),native_state('rt'))
 beads=inp('streptavidin','DynabeadsC1;25°C45min',M.SELECTION);bound=prod('bound','beadboundbiotinRNA/cDNA hybrids pluspossibleofftargets','rt');unbound=g.product('unbound','removedcapturesupernatant',M.WASTE)
 g.react('capture','CapturebiotinylatedRNA/cDNA hybrids',K.AFFINITY_CAPTURE,((x,I.TEMPLATE),(beads,I.SELECTION)),(bound,unbound));g.discard(unbound,'Physicalunboundfraction; not a guaranteeofcompleteunmodifiedRNA exclusion.')
 x=clean('capture_washes','Five4MSHAPEwashes then twoPBSwashes',bound)
 rnase=inp('rnase_release','RNaseA/T1/H plusD-biotin;37°C30min',M.CATALYST);p3=inp('p3_primer','P3primer included in sourceelutionmix; concentrationdiscrepancy unresolved',M.OLIGO,Construct([Segment('P3',P3)]))
 x=op('digest_rna','DigestoriginalRNA oncapturedbeads; not physicaldisposal',K.DIGESTION,x,((rnase,I.CATALYST),(p3,I.STRUCTURAL)),native_state('released'))
 heat=inp('dmso_heat','DMSO95°C4min; magneticbeadseparation',M.CONDITION)
 released=prod('released','elutedactualssDNA pluspossibleRNAhydrolysates','released');beadw=g.product('bead_waste','magnetically separatedbeads',M.WASTE)
 g.react('elute','DMSO/heat release and magnetictransferofsupernatant',K.STRAND_SEPARATION,((x,I.TEMPLATE),(heat,I.CONDITION)),(released,beadw));g.discard(beadw,'Sourcebeadseparation; cDNA-containing supernatant retained.')
 x=clean('cdna_column','ZymoDNA C&C5 nonstandard7×buffer+ethanol;recoverSS≥20nt',released);x=clean('cdna_page','6%denaturingPAGE isolatecDNA>70nt;SpinX/Amicon',x);x=clean('cdna_page_column','ZymoDNA C&C5 sourceSSconditions',x)
 cl=inp('circligase','CircLigaseII/Mn;60°C2h',M.CATALYST);circle=prod('circle','actual sourceinterruptedssDNA circle','circle')
 g.react('circle_closure','Circularize native5′P/3′OH cDNA',K.CIRCULARIZATION,((x,I.TEMPLATE),(cl,I.CATALYST)),(circle,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='Sourceprimer5′P andactualRTcDNA3′OH; internalchemicalspacers remainonoriginalcircle.'))
 x=clean('circle_cleanup','ZymoDNA C&C5 sourceSSconditions',circle)
 p5=inp('p5_primer','sourceP5-Solexa',M.OLIGO,Construct([Segment('P5',P5)]));pol=inp('phusion','PhusionHF+SYBR;realtimeqPCR untilsignal8000–10000;NAI-N3expected8–11cycles',M.CATALYST)
 lib=prod('pcr_library','PCRdaughters plusnativecircle/linear/offtargetclasses','library',M.LIBRARY)
 g.react('library_pcr','Amplify actualnativecircularcDNA withsourceP5/P3',K.PCR,((x,I.TEMPLATE),(p5,I.STRUCTURAL),(p3,I.STRUCTURAL),(pol,I.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(x.state,reason='Actualsourcecircle plusnativeiSp18spacers/residualclasses entersPCR; firstcycle/primerlanding/fullstrandassembly not replacedbycompletedendpoint.'))
 x=clean('pcr_column','ZymoDNA C&C5 sourcePCRconditions5×buffer;freeze-dry',lib);x=clean('library_page','6%nativePAGE cutabove~150bp;exclude120–130bpdimers;SpinX/Amicon',x);x=clean('final_column','ZymoDNA C&C5;~24µlfinal',x)
 g.sequence(x,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':75},read_uses={'Read 1':'extract_then_include'},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('cellular_rna','nai_n3','acylate','rna_isolation','isolate_rna','polya1_selection','polya1','polya2_selection','polya2','polya_rneasy_selection','polya_rneasy','dibo_biotin','click_condition','click','click_cleanup_selection','click_cleanup','fragmentation_condition','fragment','fragment_stop','stop_fragment','fragment_cleanup_selection','fragment_cleanup','end_repair','repair','adapter','t4_rnl1','ligate','ligation_cleanup_selection','ligation_cleanup','rna_page_selection','rna_page','rna_page_column_selection','rna_page_column','rt_primer','ssiii','rt','streptavidin','capture','capture_washes_selection','capture_washes','rnase_release','p3_primer','digest_rna','dmso_heat','elute','cdna_column_selection','cdna_column','cdna_page_selection','cdna_page','cdna_page_column_selection','cdna_page_column','circligase','circle_closure','circle_cleanup_selection','circle_cleanup','p5_primer','phusion','library_pcr','pcr_column_selection','pcr_column','library_page_selection','library_page','final_column_selection','final_column','final_column_out')
