"""Published D-seq WT borohydride/RT-stop route; source-owned private partial graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,Terminal,Disposition,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
SOURCE='PublishedPMC9129914 mainMethods155/159–164,Results108,actualFig1PDFp3. Oneoriginal; SIoriginals/vendor/manuals unreviewed, no confirmedabsence.'
NOTICE='PUBLISHEDSOURCE DIFFERENCE: Fig1D saysRNA100–120nt selection, mainMethods159says70–80nt. ModelselectsexplicitMethods70–80nt; exact relation between figurescale andMethodsRNAwindowunresolved. Full-length90–100ntcDNAselection belongs separateRNA-seqcontrol.'

def unknown(label):
 sc=Scene();sc.strand('unassembled',[Segment('unassembled','['+label+']',placeholder=True)]);return sc

def first_interface(*,inferred=True):
 sc=Scene();sc.strand('rna',[Segment('RNA','[selectedRNAsequenceunassembled]',placeholder=True,inferred=inferred)],mod3='OH');sc.strand('adapter',[Segment('RNA_adapter','[sourcepreadenylatedadapterbodyunassembled]',placeholder=True)],mod5='App')
 sc.ligation_end('adapter','RNA_adapter',kind='single-strand',polarity='5-prime',owner='adapter',text='Sourcepreadenylatedadapter5App')
 sc.ligation_end('rna','RNA',kind='single-strand',polarity='3-prime',owner='substrate',text='Inferredligation-compatibleRNA3OHclass')
 sc.note('rna','INFERRED ligation-compatibleRNA3OHclass, notallendprocessedRNA. Source3adapter/preactivationknown; blocked3end/wholeordering andsealedjunctionunassembled.');return sc

def protocol_graph():
 g=ProtocolGraph('D-seq selectedWTmRNA borohydrideroute',authority=SOURCE,representation='published_instance')
 x=g.input('poly_a_rna','SourceWTpoly(A)+mRNAaftertwo oligo-dTcellulosecaptures; upstreamhotacidphenol/extraction/8mgtotalRNAinputnotassembled. No>200ntLCMS-specificcolumnimported.',K.MOLECULAR,primary=True)
 def op(key,label,kind,items,product,geometry=None,reason=None,extra=()):
  nonlocal x
  ins=[(x,R.TEMPLATE)]
  for ident,text,mk,role in items:ins.append((g.input(ident,text,mk),role))
  y=g.product(key+'_product',product,K.MOLECULAR);kw={}
  if geometry is not None:kw['geometry']=geometry
  if reason:kw['reason']=reason
  g.react(key,label,kind,ins+list(extra),(y,),**kw);x=y
 def c(key,text):return(key,text,K.CONDITION,R.CONDITION)
 def e(key,text):return(key,text,K.CATALYST,R.CATALYST)
 def clean(key,label):
  nonlocal x
  b=g.input(key+'_condition',label,K.CONDITION);kept=g.product(key+'_kept','Sourcephysicallyretainedpreparation;actualrecovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Physicalprecipitation/cleanupfractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(kept,w));g.discard(w,'Physicalworkup,notcomputationalfilterorperfectselection.');x=kept
 def size_select(key,label):
  nonlocal x
  gel=g.input(key+'_gel',label,K.SELECTION);kept=g.product(key+'_kept','Sourcegel-selectedmaterial; exactrecoveryandallcompanionstatesunassembled',K.MOLECULAR);w=g.product(key+'_outside','Physicalothergelregions/workupfractions',K.WASTE);g.react(key,label,S.SELECTION,((x,R.TEMPLATE),(gel,R.SELECTION)),(kept,w));g.discard(w,'Sourcephysicalsizefractionation; notsequencereadfiltering.');x=kept
 op('fragment','Source10mMZnCl2,94C1min',S.FRAGMENTATION,(c('zinc_condition','Sourcezincfragmentationconditions;initialendsunassembled'),),'SourcefragmentedRNAandchemicalcompanions;exactenddistributionunassembled');clean('fragment_precip','Sourceprecipitationafterzincfragmentation')
 op('reduce','Source18ulRNAwater+2ul10mgmlNaBH4in500mMTrispH7.5stock,0C1h',S.CHEMICAL_CONVERSION,(c('bh4','Sourceborohydride/Trisstock;500mMnotfinal20ulbufferconcentration'),),'SourceD→tetrahydrouridine-treatedRNA,unconvertedDandotherreactivebases/sideproducts; noabsoluteDselectivityorcompleteconversion')
 op('quench','Source4ul6Naceticacid',S.QUENCHING,(c('acid','SourceCH3COOHquench'),),'SourcequenchedRNAreaction; degradationproductsnotphysicallygone');clean('reduce_precip','Sourceprecipitationafterborohydridequench')
 op('end_process','SourceRNAdephosphorylation/endrepairwithCIPandPNK;actualadditionorder/ATP/wholeformulationunassembled',S.MAGIC,(e('cip','SourceCIP'),e('pnk','SourcePNK')),'SourceprocessedRNAendclasses;unprocessedandfailedclassesremain',reason='SourceCIP/PNKpositive, exactadditionprogramme/endpopulationsnotassembled;nosilentinferenceallRNA5P/3OH.')
 size_select('rna_select','SelectedMethods70–80ntRNA8percentureaTBEgel;Fig1D100–120ntdiffers. Overnight4C400ulelution300mMNaOAc5.5/1mMEDTA/100UmlRNasin')
 clean('rna_precip','SourceprecipitationafterRNAgelelution')
 a=g.input('three_adapter','SourcepreadenylatedIDTRNA3adapter;wholebody/modificationsunassembled',K.ADAPTER,chemical_state={'source_activation':'preadenylated','terminal_blocking':'unassigned'})
 geom=AdapterLigationGeometry.known('single-strand',first_interface(),statement='SourceRNA3adapterligation/preactivationknown;ligatableRNA3OHclassisinferredandpartial. Exactadapterbody/blocking/stateandsealedjunctionunassembled.')
 op('rna_ligate','SourceT4RNAligaseinATP-freebuffer50mMTris7.8/10mMMg/10mMDTT,22C2.5h',S.LIGATION,(e('rna_ligase','SourcefirstT4RNAligaseNEB;isoformunassigned,notborrowedfromsecondreaction'),c('rna_lig_buffer','SourceATP-freebuffer')),'RNA3adapterligationmixture,unligatedRNA/adapterandbyproductsretained',geometry=geom,extra=((a,R.STRUCTURAL),))
 clean('rna_ligate_precip','SourceprecipitationafterRNAadapterligation')
 primer=g.input('rt_primer','SourceRTprimer,wholeordering/modificationsunassembled',K.OLIGO);b=g.input('rt_anneal_condition','Source95C5minthenice5minRNA/primeranneal',K.CONDITION);template=g.product('rt_template','Same-tubeRNAtemplatecomponentafterannealing',K.MOLECULAR);ret=g.product('rt_retained_primer','Co-retainedexternalRTprimercomponent;notanotherstockaddition',K.OLIGO);g.react('rt_anneal','SourceRNA/primerannealing',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(b,R.CONDITION)),(template,ret))
 rt=g.input('ssiii','Source1ulSSIII',K.CATALYST);nts=g.input('rt_dntps','Source1ul10mMdNTPsubstrates',K.MOLECULAR);buf=g.input('rt_conditions','Source1ul.1MDTT/1ulRNasin/3ul5xFirstStrandBuffer',K.CONDITION);y=g.product('rt_product','SourcecDNA/RNA/primerreactionincludingD-associatedstopsandfull-lengthproducts;notallstopsassignedtoD',K.MOLECULAR);g.react('rt','Source50C1hSSIII',S.REVERSE_TRANSCRIPTION,((template,R.TEMPLATE),(ret,R.STRUCTURAL),(nts,R.STRUCTURAL),(rt,R.CATALYST),(buf,R.CONDITION)),(y,));x=y
 size_select('cdna_select','Source50–80nttruncatedcDNA8percentureaTBEgel;400ul300mMNaCl/10mMTris8overnightRTelution;elution/precipitationcompoundworkupasstated')
 used=g.product('cdna_5ul','Source5ulelutedcDNAselectedforlinkerligation',K.MOLECULAR);rest=g.product('cdna_remainder','UntransferredDNAeluate;actualrecoveredvolumeunassigned',K.MOLECULAR);g.react('cdna_aliquot','Source5ulactualtransferfromgel-elutedstock',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No400ulperfectrecovery/395ulremainderassumed.'));x=used
 linker=g.input('cdna_linker','Source0.8ul80uMIDTlinkernamed5adapter;actualchemicalend/wholeorderingunassembled',K.ADAPTER);mix=g.input('linker_anneal_condition','Source1ulDMSO/75C2minthenice2min',K.CONDITION);bound=g.product('linker_preparation','Same-tubecDNAcomponentwithco-retainedlinker;fullgeometryunassembled',K.MOLECULAR);link=g.product('retained_linker','Co-retainedIDTlinkercomponent;noextraaddition',K.ADAPTER);g.react('linker_prepare','SourcecDNA/linkerDMSOthermalpreparation',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(linker,R.STRUCTURAL),(mix,R.CONDITION)),(bound,link));x=bound
 lig=g.input('cdna_ligase','Source0.5ulT4RNAligase2;variantunassigned',K.CATALYST);buff=g.input('cdna_lig_condition','Source2ulligasebuffer/.2ul.1MATP/6.5ulPEG8000/3.6ulwater,overnight22C;statedvolumesnotnormalized',K.CONDITION);y=g.product('cdna_ligate_product','SourcecDNA/linkerligationmixtureandunligated/off-productclasses;exactpolarity/sealedproductunassembled',K.MOLECULAR)
 g.react('cdna_ligate','Sourcecalled5adapterligationonRT-stopcDNA;notautomaticallyassignedcDNA5end',S.LIGATION,((x,R.TEMPLATE),(link,R.STRUCTURAL),(lig,R.CATALYST),(buff,R.CONDITION)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown('SourcecDNA/linkerligationendsunassembled'),),'Actualsource5adapterterm/ligase2/ATPmixerretained;whichchemicalendandwholejunctiongeometrynotassembled,nosourceabsenceclaim.'));x=y
 clean('silane','SourceDynabeadMyOneSilane cleanupmanufacturerboundary')
 p1=g.input('pcr_forward','SourcePCRforwardprimerwholeorderingunassembled',K.OLIGO);p2=g.input('pcr_reverse','SourcePCRreverseprimerwholeorderingunassembled',K.OLIGO);pol=g.input('phusion','SourcePhusionDNApolymerase',K.CATALYST);y=g.product('pcr_product','Sourceampliconandoriginaltemplate/freeprimer/backgroundmixture;wholeproductunassembled',K.MOLECULAR);g.react('pcr','SourcePhusionPCRfrombothlinkers;cycles/unassignedsampleindexesnotinvented',S.PCR,((x,R.TEMPLATE),(p1,R.STRUCTURAL),(p2,R.STRUCTURAL),(pol,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(unknown('SourcewholePCRtemplate/primersunassembled'),reason='Wholeoperationalprimeridentitiesandlandingunassembled;genotypingprimersnotsubstituted.'));x=y
 size_select('pcr_gel','SourcePCRproductgelpurification;exactproductionwindowunassigned');clean('pcr_precip','SourceprecipitationafterPCRgelpurification')
 other=g.input('other_libraries','ActualotherPCRproductsinthestudypool;number/compositionunassembled',K.MOLECULAR);pool=g.product('pool','Sourcepooledsequencingpreparation',K.MOLECULAR);g.react('pooling','SourcephysicalPCRlibrarypool',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 lib=g.product('library','SourceD-seqHiSeq2500librarypool;wholeindex/primer/readgeometryunassembled',K.LIBRARY);g.react('loading','Physicallibraryidentity',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceIlluminaHiSeq2500;runlength/indexcycles/wholeprimersunassigned. RTstops1nt3primeofDareobservedconditionalreadout,notallRNAcleavedatD. '+NOTICE));g.add_source_item(SourceItem('pending','Wholeadapter/RT/PCRends,junctions,purificationandreadlayout','Selectedpublishedmainonly;remainingSI/linkedoriginalsunreviewed,noconfirmedabsence.'));return g

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph()
 c('NoTAPorDUSgeneticconstructioninvented',not any('tap'in s or'dus'in s for s in g.steps),True)
 c('RNAendrepairkeepsbothsourceenzymes',[i.material for i in g.steps['end_process'].inputs if i.role==R.CATALYST],['cip','pnk'])
 c('Firstligationconsumesactualsourceadapter','three_adapter'in[i.material for i in g.steps['rna_ligate'].inputs if i.role==R.STRUCTURAL],True)
 c('RTconsumesproducedprimer',[i.material for i in g.steps['rt'].inputs if i.role==R.STRUCTURAL and g.materials[i.material].kind==K.OLIGO],['rt_retained_primer'])
 c('Actual5ulstockremainderretained','cdna_remainder'in g.terminals,True)
 c('Secondligationusesretainedlinker',[i.material for i in g.steps['cdna_ligate'].inputs if i.role==R.STRUCTURAL],['retained_linker'])
 c('Knownpartialfirstendandunassembledsecondend',(g.steps['rna_ligate'].geometry.availability,g.steps['cdna_ligate'].geometry.availability),('known','unavailable'))
 src=Source(sources_root()/'d-seq-draycott2022__10.1371+journal.pbio.3001622'/'journal.pbio.3001622_PMC9129914.xml.txt','PublishedmainMethods159/Fig1Doriginalp3')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('SourceRNAselectionMethods70to80','RNAfragments(70to80nt)'in t,True);c('SourcecDNAtruncationwindow50to80','TruncatedcDNAsweresizeselected(50to80nt)'in t,True)
 c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
