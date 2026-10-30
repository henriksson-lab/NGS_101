from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s
AUTHORITY='Published2016PMC4889966 mainMethods120/Results36–50; actualFig1PDFp2/Methodsp7. ONEoriginal; inaccessibleactualDC1landing/citedIossifov2012ordinaryunread, noSMASH2recipeinheritance.'
def qc(b,key,label):
    g=b.g;y=g.product(key+'_analytical',label,K.MOLECULAR);r=g.product(key+'_stock','Actualsourcepreparationremainingafteranalyticalallocation; amountunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(y.id,Disposition.QC,label));b.x=r

def keep_supernatant(b,key,label,beadlabel):
    g=b.g;y=g.product(key+'_supernatant','Actualretainedsupernatant; recovery/sizeselectivityunassigned',K.MOLECULAR);w=g.product(key+'_bead_fraction','Actualremovedbead-boundfraction',K.WASTE);g.react(key,label,S.CLEANUP,((b.x,R.TEMPLATE),b.stock(key+'_beads',beadlabel,K.SELECTION,R.SELECTION),b.cond(key+'_condition','SourceRT5min; physicalfractionselectionnotperfectsizepurity')),(y,w));g.discard(w,'Actualbead-boundfractionexcludedbyselectedsourceoperation');b.x=y

def polish(b,key):
    b.op(key,'SourceendpolishRT30min',S.MAGIC,(b.stock(key+'_t4pol','ActualT4DNApolymerase'),b.stock(key+'_klenow','ActualDNApolymeraseIlargeKlenowfragment'),b.stock(key+'_pnk','ActualT4PNK'),b.stock(key+'_supplies','Actualendpolishbuffer/nucleotide/cofactorsupplies; identities/doses/packagingunassigned',K.MOLECULAR,R.STRUCTURAL)),reason='Sourcepositiveendpolish; actualheterogeneousfragmentends/filledproducts and reagentformulationsunassembled. No inferredperfectblunt5P/3OHpool.')

def ligase_supply(b,key):
    g=b.g;st=g.input(key+'_kit','ActualQuickLigationKitallocation; formulation/doseunassembled',K.MOLECULAR);enzyme=g.product(key+'_ligase','Ligationcatalystcomponentofsameactualkitallocation; completeformunassigned',K.CATALYST);other=g.product(key+'_supplies','Otheractualsamekitbuffer/cofactorcomponents; identities/dosesunassigned',K.MOLECULAR);g.react(key+'_parts','Copresentcomponentsofoneactualkitallocation; no extraexternalenzyme',S.GENERIC,((st,R.TEMPLATE),),(enzyme,other),changes_molecule=False);return((enzyme,R.CATALYST),(other,R.STRUCTURAL))
def clean(b,key,label):b.clean(key,label,extra=(b.stock(key+'_selection','Actual'+label+'selectionmaterial',K.SELECTION,R.SELECTION),))
def ligation_geometry():
    sc=Scene();sc.strand('substrate',[Segment('Unassembledrepresentativegenomicfragment','[unassembledfragment]',placeholder=True,inferred=True),Segment('Source3primeAclass','A')]);sc.ligation_end('substrate','Source3primeAclass',kind='cohesive',polarity='3-prime',owner='substrate',text='Sourcepositiveintended3A-fragmentclass; wholefragment/endsdistributionunassembled, notperfectAtailing')
    ad=Scene();ad.strand('adapter',[Segment('Wholecitedadapterunassembled','[unassembledadapter]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,ad),'SourcepositiveA-tailedsubstrateclass/local3A retained; citedIossifovadapter/SIunread, wholeadapter/endmods/occupiedjunctionsunassembled. Noabsenceclaimorcanonicaladapterinherited.')

def protocol_graph():
    b=B('SMASHinitialSSC-familysonication/restrictionpartialpreparation',AUTHORITY,'Actualselected1ugwholeblood-derivedSSCgenomicDNAin75ul10mMTris8.5');g=b.g
    b.op('shear','SourceCovarisE210shearingtoabout100bp',S.MAGIC,(b.cond('covaris','ActualsourceCovarisprogramme delegatedmanufacturerparametersunassembled'),))
    b.op('restriction','CviKI-1+NlaIII37C1hin90ulCutSmart',S.DIGESTION,(b.stock('cviki','ActualCviKI-1'),b.stock('nlaiii','ActualNlaIII'),b.cond('restriction_buffer','Actual1xCutSmart90ulreaction; enzymeallocationsunassigned')))
    b.op('concentrate','SpeedVacreducesvolumetoabout30ul',S.GENERIC,product='ActualconcentrateddigestedDNAstock; volumerecoveryunassigned')
    keep_supernatant(b,'short_select','Source2.5xAMPure removeslongfragments; RETAINsupernatant','Actual2.5xAMPureallocation')
    clean(b,'short_column','QIAquicknucleotideremovalkit/elute30ulwater');qc(b,'short_qc','SourceBioanalyzerfragmentmeasurement; source40–50bpaverage, analyticalamountunassigned')
    polish(b,'polish1');clean(b,'polish1_clean','QIAquicknucleotideremovalkit/elute30ulwater')
    b.op('concatenate','SourcerandomintermolecularligationRT15min',S.MAGIC,ligase_supply(b,'concat'),product='ActualchimericDNApreparationwithphysicallyjoinedindependentgenomicintervals pluslinear/unligated/offproductclasses; wholejunctionsunassembled',reason='Positivephysicalfragmentligation, notchromatinproximityorcomputationalconcatenation. Wholepolishedends/fragmentidentities/circularandoffproductfatesunassembled.')
    clean(b,'concat_clean','1.6xAMPureXP');polish(b,'polish2')
    b.op('a_tail','SourceKlenowAaddition37C30min',S.MAGIC,(b.stock('tail_klenow','PrintedKlenowfragment(3prime→5primeexo); sourceenzymeformunassigned'),b.stock('tail_supplies','ActualA-tailingreactionbuffer/nucleotidesupplies; identities/amounts/packagingunassigned',K.MOLECULAR,R.STRUCTURAL)),reason='SourcepositiveadditionofsingleAto3ends; actualPDFp7printsKlenowexo withoutminus. No inferredexo-minusvendorrecipe orcompleteenddistribution.')
    clean(b,'tail_clean','1.6xAMPureXP')
    barcode=feature('smash.sample','sample_index','unknown',note='Sourcebarcodedsequencingadapters/samplemultiplexing; fullorders/orientation/widthunassembled, notUMI')
    ad=Scene();ad.strand('source',[Segment('Sourcebarcodedadapter','[wholebarcodedadapterunassembled]',placeholder=True,feature=barcode)])
    b.op('adapter_ligate','SourcebarcodedadapterligationbyQuickLigation',S.MAGIC,ligase_supply(b,'adapter_kit')+(b.stock('adapter','ActualIossifov2012-delegatedbarcodedadapterpreparation; constituentwholeordersunassembled',K.ADAPTER,R.STRUCTURAL,state=ad),),geometry=ligation_geometry())
    clean(b,'adapter_clean','1.6xAMPureXP/elute50ulwater')
    keep_supernatant(b,'large_reject','Source30ulbeads(.6x)into50ulDNA; RETAINsupernatant','Actual30ulAMPureXP')
    clean(b,'target_capture','Additional8ulAMPure(0.16xORIGINAL50ul),retainbead-boundfraction/two180ul80percentEtOH/elute30ulwater')
    b.split('pcr_select','ActualselectedDNAallocationfor8cyclePCR; amountunassigned','ActualPCRtemplateallocation','Untransferredsize-selectedstock')
    stock=g.input('phusion_stock','Actual1xPhusionHighFidelityPCRMasterMixwithHFbuffer allocation; volumeunassigned',K.MOLECULAR);pol=g.product('phusion_enzyme','PolymerasecomponentofsameactualPhusionmastermix',K.CATALYST);nt=g.product('phusion_nt','Nucleotidecomponentsofsameactualmastermix',K.MOLECULAR);buf=g.product('phusion_buffer','HFbuffercomponentofsameactualmastermix',K.CONDITION);g.react('phusion_parts','Copresentcomponentsofoneactualmastermixallocation; no freshindependentadditions',S.GENERIC,((stock,R.TEMPLATE),),(pol,nt,buf),changes_molecule=False)
    primers=g.input('pcr_primers','ActualPCRprimer supplies aggregateinventory; exactbodies/speciescount/packaging/stocknumbersunassigned',K.MOLECULAR)
    f=g.product('pcr_forward','ActualfunctionalforwardPCRprimerclass; wholeorders/speciescountunassembled',K.OLIGO,state=unknown_scene());r=g.product('pcr_reverse','ActualfunctionalreversePCRprimerclass; wholeorders/speciescountunassembled',K.OLIGO,state=unknown_scene());g.react('pcr_primer_parts','Exposefunctionalcomponents ofactualaggregatePCRprimerinventory; not claimedcopackagingoradditionalfreshdoses',S.GENERIC,((primers,R.TEMPLATE),),(f,r),changes_molecule=False)
    b.op('pcr','Sourceeight-cyclePCR',S.PCR,((pol,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Sourceactual8PCRcycles; fullprimers/amplicons/occupiedsitesandthermocyclingunassembled. Functionalforward/reverseclasses doNOTassignspeciescounts/packagingorcanonicalsequence.'))
    qc(b,'library_qc','ActualBioanalyzerquantificationanalyticalallocation; amountunassigned');b.op('normalize','Sourcedilutionto10nM',S.GENERIC)
    b.split('sample_pool_select','Actualsampleallocationforindexedsequencingpool; amountunassigned','Selectedsampleshare','Untransferrednormalizedsamplestock')
    b.op('pool','Sourceindexedmultiplexingpool; otheractualsampleallocationsunassembled',S.MAGIC,(b.stock('other_samples','Otheractualbarcodedlibraryallocations; count/volumes/sourceconstructsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.split('loading','Actualloadingaliquotfrompreparedpool; amountunassigned','Selectedpooledloadingaliquot','Untransferredpreparedsequencingpool')
    lib=g.product('library','SourceSSCinitialSMASHsequencinglibrary; fullchimeras/adaptergeometryunassembled',K.LIBRARY);g.react('load','SelectedSSCHiSeq2000PE101; SKBR3PE150outsidechosenroute',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':101,'R2':101},limitation='PublishedinitialSSCbranch; wholeadapter/index/PCRprimerorders/landingunassembled. SMASHmapparsing/duplicatefiltering iscomputational, notDNAcleavage.'))
    g.add_source_item(SourceItem('size_scope','PublishedFig1caption250–700bp; Results/Methods300–700bp.','Actualbeadfractionoperationsrepresentedwithoutclaimingperfectcutoffs; no silent250/300repair.'))
    g.add_source_item(SourceItem('source_scope',AUTHORITY,'Privatepartialsourcepositivepreparation; no predecessoradaptersequenceinheritance/fullchemistry/acceptedabsence.'));return g
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('EarlyshortfragmentselectionRETURNSsupernatant',g.steps['short_column'].inputs[0].material,'short_select_supernatant')
    c('Finalsecondbeadadditionusesfirstsupernatant',g.steps['target_capture'].inputs[0].material,'large_reject_supernatant')
    c('Physicalchimeraproductionbeforeadapterligation',g.steps['concat_clean'].inputs[0].material,'concatenate_product')
    c('Secondpolishusescleanedchimerastock',g.steps['polish2'].inputs[0].material,'concat_clean_kept')
    c('Phusionenzymefromactualsamepremix',g.producer_of('phusion_enzyme').id,'phusion_parts')
    c('PCRusesselectedfractionnotwholepreparedstock',g.steps['pcr'].inputs[0].material,'pcr_select_selected')
    c('AnalysisQCexcludedfromproduction',g.terminals['library_qc_analytical'].disposition,Disposition.QC)
    c('Onlyloadingaliquotsequenced',g.steps['load'].inputs[0].material,'loading_selected')
    c('Wholeadapterunassembled',g.materials['adapter'].state.strands['source'].segs[0].placeholder,True)
    c('PCRforwardandreversecomefromoneaggregateinventory',g.producer_of('pcr_forward').id==g.producer_of('pcr_reverse').id=='pcr_primer_parts',True)
    c('PCRstepistypedwithunassembledgeometry',g.steps['pcr'].kind,S.PCR)
    c('KnownAtailinterfacebelongs tosubstratenotunknownadapter',[(e.owner,e.polarity)for sc in ligation_geometry().scenes for e in sc.ligation_ends],[('substrate','3-prime')])
    c('SourceSSCrunselected',g.sequencing['library'].read_lengths,{'R1':101,'R2':101})
    c('Privatepartialnotpublicationready',bool(g.publication_problems()),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'gr.201491.115_PMC4889966.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1101/gr.201491.115')):
        t=p.read_text()
        for name,needle in [('Originalsecondadditiondenominator','0.16× the original 50 µL'),('Originalphysicalconcatligation','randomly ligated to form longer stretches of chimeric DNA'),('OriginalPhusionHFmix','Phusion high-fidelity PCR master mix with HF buffer'),('OriginalSSCrun','HiSeq 2000 (2 × 101 bp run mode')]:c(name,needle in t,True)
    c.report()
if __name__=='__main__':selftest()
