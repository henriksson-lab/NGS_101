"""Original2018 SCOPE-seq; source-ownedfactory/DNAlookup/cellularRNA partial graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,K.MOLECULAR,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source copresent reaction mixture; incomplete/unreacted classes and carryover unassigned'):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Positive source operation; historical component formulation/whole molecular products unassembled, not source absence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def split(self,key,label,branches):
        yy=tuple(self.g.product(key+'_'+name,txt,K.MOLECULAR)for name,txt,disp in branches);self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),yy)
        for y,(_,txt,disp)in zip(yy,branches):
            if disp is not None:self.g.mark_terminal(Terminal(y.id,disp,txt))
        self.x=yy[0];return yy
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Source retained physical fraction; recovery/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated wash/workup fraction',K.WASTE);self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source separated workup; no perfect specificity/yield.');self.x=y
from dataclasses import dataclass
from chemdraw import feature
from illumina import P5,P7
from nextera import ADAPTOR_S7
from pcr import PCRGeometry
CODES=('GAGTGATTGAACGCCCTC','TGGGTGCTTCGTTATTCC','CTGTGTAGAGGGCGAGAC','GTTTGAAGCACACACACC','CTCGCAGGTACAACATCA','CTGGAATGTCCCTGGTTC','CGTGATCCCGTAAATGCT','TGGCAGACATTCCCAATC','TTGTGGAGTTGAGGAGGG','CGAGAGTCAAATCAGGCA','TTTCACAGATCTCAGCGG','TACTGAAGAAGCCTTGCC')
OBO_PREFIX='TTTTTTUTTTT'
READ1='GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC'
UNIVERSAL=P5+READ1
UNIVERSAL_RAW=UNIVERSAL[:-2]+'*'+UNIVERSAL[-2]+'*'+UNIVERSAL[-1]
NEXTERAPOLYA=ADAPTOR_S7+'A'*30
SMRT='AAGCAGTGGTATCAACGCAGAGT'
LOOKUP_P7=P7+'CGTGTTCA'+'GTCTCGTGGGCTCGG'
@dataclass(frozen=True)
class OpticalOligo:
    code:str
    prefix:str=OBO_PREFIX
    spacer:str='/iSp18/'
    terminal:str='C'
    @property
    def ordering(self):return '/5Phos/'+self.prefix+self.spacer+self.code+'/3ddC/'
    @property
    def nucleotide_count(self):return len(self.prefix)+len(self.code)+len(self.terminal)
    @property
    def parts(self):return (('DNA',self.prefix),('non_nucleotide',self.spacer),('DNA',self.code),('modified_DNA',self.terminal,'3ddC'))
    def __post_init__(self):
        if len(self.code)!=18 or set(self.code)-set('ACGT'):raise ValueError('Actualsource18nt opticalcode required')
        if self.spacer!='/iSp18/'or self.terminal!='C':raise ValueError('Thislocaldescriptoronlycoversliteraloriginalsourcechemistry')
def oligo_scene(seq,mod5='',mod3=''):
    s=Scene();s.strand('ordering',[Segment('whole_ordering',seq)],mod5=mod5,mod3=mod3);return s
def lookup_p7_scene():
    s=Scene();s.strand('primer',[Segment('p7',P7),Segment('index','CGTGTTCA',feature=feature('scope.lookup.index1','sample_index','fixed')),Segment('annealing','GTCTCGTGGGCTCGG')]);return s
def optical_pair(n):
    s=Scene();seg='code'+str(n+1);s.strand('code',[Segment(seg,CODES[n],feature=feature('scope.optical.'+str(n+1),'feature_barcode','fixed',group='scope.optical',note='Onecomponentofbead opticalcombination, notRNAUMI'))],label='PartiallocalOBOcodeview; upstreamSp18/DNAprefixanddownstreamddC retainedinchemical_state')
    s.anneal('probe',[Segment('probe',revcomp(CODES[n]))],to='code',pair=('probe',seg),mod5='Cy5');return s
def unknown_scene():
    s=Scene();s.strand('unknown',[Segment('unassembled','[whole source product unassembled]',placeholder=True)]);return s
def protocol_graph():
    b=B('Original SCOPE-seq2018 beadfactory/lookup/RNA','PublishedPMC6305572main39–61; publishedAdditionalFile1Methods3–8/97/99; publishedAdditionalFile4TableS2. Threecombinedoriginals, noSCOPE2/vendor/predecessorconstructorinheritance.','Sourceapproximately700,000commercialChemGenesMacosko-2011-10(V+)beads; actualcaptureoligowholebody/UMIlengthunassembled');g=b.g
    # Microplate is a batch of separate compartments, never a one-pot reaction or prematurely pooled DNA.
    lookup=[]
    for rd in [1,2]:
        beads=g.product('factory'+str(rd)+'_beads','Retainedbeadfractions in64separatewellcompartments; underlyingcapturesequencesunassembled',K.MOLECULAR)
        free=g.product('factory'+str(rd)+'_supernatants','64distinctbead-freeDNAsupernatantstocks; NOTphysicallypooled',K.MOLECULAR)
        g.react('release'+str(rd),'Source64wells/5minsonication/4C12h; transferdetachedDNAto separateplate−20C, washretainedbeads',S.MAGIC,((b.x,R.TEMPLATE),b.cond('factory_conditions'+str(rd),'Sourcewater/64wellcompartments; exactdetachmentmechanismunassembled')),(beads,free),reason='Sourcephysicalbead/supernatantseparationpositive; fullcommercialbeadchemistry andreleasebondstatesunassembled.');lookup.append(free);b.x=beads
        obos=g.input('obo_round'+str(rd),'Sourceper-wellselectedcombination4uMeachfromsixOBOs; actualwellcombinationunassigned, NOTall-sixonepotaddition',K.ADAPTER,chemical_state={'source_inventory':tuple(OpticalOligo(c)for c in CODES[(rd-1)*6:rd*6]),'well_assignments':'unassembled','stock_semantics':'64separatewell-specificaliquots'})
        b.op('ligate_optical'+str(rd),'Source64separatereactionsRnl1/.7Uul/ATP1mM/PEG25percent/Tween.1percent24hRT',S.MAGIC,((obos,R.STRUCTURAL),b.stock('rnl'+str(rd),'SourceT4RNAligase1per-wellpreparation'),b.cond('ligbuffer'+str(rd),'Source1xRnlbuffer/1mMATP/25percentPEG/.1percentTween')),product='64separatewellbead-boundligated/unligatedOBO/capture-siteclasses; source3ddC blocksfurtherextension, exactoccupancyunassigned')
        b.op('stop_factory'+str(rd),'Source5ul10percentSDS+1ul500mMEDTA5min',S.QUENCHING,(b.cond('stop'+str(rd),'SourceSDS/EDTA'),))
        b.op('pool_factory'+str(rd),'Sourcepoolbeadsfrom64wells; separatewellproductionoutlinedabove, noDNAlookupstockpooling',S.MAGIC,(b.cond('pool_conditions'+str(rd),'Physicalpoolof64retainedbeadfractions'),),product='Sourceactuallypooleddual/partialopticallytaggedbeadstock')
        b.clean('factory_wash'+str(rd),'SourceTE/TWpoolwash; opticalsiteoccupancy/stocklossunassigned')
    final_beads=b.x
    # One explicit lookup well; other127remainseparate until theirPCRaliquotpool.
    b.x=lookup[0];first,other63=b.split('lookup_well','Sourceoneof64first-roundsupernatantwells selected; identitiesunassigned',[('one','Selectedfirst-roundwell DNAstock; illustrativeP7primer1choice',None),('other63','Other63first-roundsupernatantstocks; separatecompartments',None)])
    b.x=first;b.split('lookup5','Source5ulsupernatantselected',[('selected','Source5ulDNA',None),('rest','Untransferredsupernatantstock',Disposition.UNRESOLVED)])
    np=g.input('nexterapolya','Source10ul1uMNexterapolyA',K.OLIGO,state=oligo_scene(NEXTERAPOLYA))
    b.op('lookup_anneal','Source2ul10xNEBuffer2/95C2min/coolRT',S.HYBRIDIZATION,((np,R.STRUCTURAL),b.cond('lookup_buffer','SourceNEBuffer2/thermalprogramme')))
    template=g.product('lookup_template','Same-tubebead-freeDNAtemplatecomponent',K.MOLECULAR);ret=g.product('lookup_retained_primer','PreviouslyaddedNexterapolyAprimercomponent',K.OLIGO,state=oligo_scene(NEXTERAPOLYA));g.react('lookup_components','Same-tubecomponents,nophysicalsplitornewprimeraddition',S.GENERIC,((b.x,R.TEMPLATE),),(template,ret),changes_molecule=False);b.x=template
    b.op('lookup_extend','Source2ul10mMdNTP+1ulKlenow;RT15min37C1h75C20min',S.EXTENSION,((ret,R.STRUCTURAL),b.stock('lookup_klenow','SourceKlenowFragmentM0212S'),b.stock('lookup_nt','Source2ul10mMdNTPmix; conventionunassigned',K.MOLECULAR,R.STRUCTURAL)))
    b.split('lookup4','Source4ulextensionreactionselectedforPCR',[('selected','Source4ulPCRtemplate',None),('rest','Untransferredextensionreactionstock',Disposition.UNRESOLVED)])
    p7=g.input('lookup_p7','SourceTableS2P7primer1; selectedexemplar/actualwellpositionunassigned',K.OLIGO,state=lookup_p7_scene());p5=g.input('lookup_p5','SourceuniversalP5 '+UNIVERSAL_RAW,K.OLIGO,state=oligo_scene(UNIVERSAL,mod3='PSatfinal2internucleotidebonds'))
    stock=g.input('lookup_pcr_stock','Source8ulKAPAreadyPCRmix; componentamountsunassembled',K.MOLECULAR);cat=g.product('lookup_pcr_cat','KAPAcatalystcomponent',K.CATALYST);nt=g.product('lookup_pcr_nt','NucleotidesofsamePCRmix',K.MOLECULAR);g.react('lookup_pcr_components','OnepreparedsourcePCRmix, notnewNTaddition',S.GENERIC,((stock,R.TEMPLATE),),(cat,nt),changes_molecule=False)
    b.op('lookup_pcr','Source2ul2uMeachP7/P5+8ulKAPA;14cycles95/55/72',S.PCR,((p7,R.STRUCTURAL),(p5,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Fullsourcewholeproducts/primerlandingunassembled, notpublicationabsence.'))
    b.split('lookup_pcr5','Source5ulperPCRforpool',[('selected','Source5ulPCRproductforpool',None),('rest','UntransferredselectedwellPCRstock',Disposition.UNRESOLVED)])
    selected_lookup=b.x
    otherpool=g.product('other_lookup127_aliquots','127separatelyindexedlookupPCR5ulaliquots; atthisoutputstillseparate',K.MOLECULAR);otherrest=g.product('other_lookup127_remainders','UntransferredDNA/extension/PCRstocksforother127separatewells',K.MOLECULAR)
    g.react('other_lookup127','Actualother127wells:5ulDNA/10ul1uMpolyA/Klenow;4ulPCR/14cycles;5ulaliquot, fullperwellstatesunassembled',S.MAGIC,((other63,R.TEMPLATE),(lookup[1],R.TEMPLATE),b.stock('other_lookup_preparation','Sourceother127actualseparateKlenow/PCRpreparedreagents/primers; wellassignmentsunassembled',K.MOLECULAR,R.STRUCTURAL)),(otherpool,otherrest),reason='Sourceparallelpreparationpositive; do notpoolDNAtemplatesorimport128indexchoices.');g.mark_terminal(Terminal(otherrest.id,Disposition.UNRESOLVED,'Actualuntransferredparallelstocks.'))
    lp=g.product('lookup_pool','Actual128lookupPCRaliquotpool',K.MOLECULAR);g.react('lookup_pool128','Sourcepool5ulperreaction AFTERindexedPCR',S.POOL,((selected_lookup,R.TEMPLATE),(otherpool,R.TEMPLATE)),(lp,));b.x=lp
    b.clean('lookup_spri','SourceSPRI2:1')
    b.op('lookup_gel_run','Source10percentTBEpolyacrylamidegelelectrophoresis BEFOREstaining',S.GENERIC,(b.cond('lookup_gel','Source10percentTBEgel/electrophoresis'),),product='SourceelectrophoreticallyseparatedDNAclassesinactualgel; desired bandnotyetcut')
    b.op('lookup_stain','Source1xSYBRGoldDNAgelstain',S.GENERIC,(b.stock('sybr','SourceSYBRGoldphysicaldye',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('lookup_page','SourceafterSYBRGoldstainingextract~180bpband/gelpurification; unselectedgelphysicalfractions')
    b.split('lookup_qc','SourceBioanalyzerQCaliquot vsstock; volumesunassigned',[('library','LookupDNAproductionstock',None),('qc','SourceBioanalyzerQCaliquot',Disposition.QC)])
    b.split('lookup_load','Sourcequalifiedlookupaliquotforloading',[('aliquot','SelectedlookupDNAaliquot',None),('rest','UntransferredlookupDNAstock',Disposition.UNRESOLVED)])
    phix=g.input('lookup_phix','Source40percentPhiXDNAcontrollibrary contribution',K.MOLECULAR);lib=g.product('lookup_library','SourceDNAlookup+40percentPhiXloadingmixture',K.LIBRARY);g.react('lookup_phix_mix','Sourceactual40percentPhiXspike',S.GENERIC,((b.x,R.TEMPLATE),(phix,R.SPIKE_IN)),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceNextSeq50026baseR1/8baseI1 customRead1 '+READ1+'; barcode12nt/readlandingwholeproductsunassembled. LookupisbeadbatchtechnicalDNA, notsecondbiologicalmodality.'))
    # CellularRNA branch consumes anactualaliquotofretainedfactorybeads.
    b.x=final_beads;production,rest=b.split('bead_aliquot','Source500beads/ulloadingaliquot; actualvolumeunassigned',[('production','SelectedfactorybeadsforcellularRNAcapture',None),('rest','Untransferredpreparedbeadbatchstock; storage/lateruseunassigned',Disposition.UNRESOLVED)])
    cells=g.input('cells','SourceU87/3T3singlecellsforliveimage/microwellRNAcapture; actualmixproportionsunassigned',K.BIOLOGICAL,primary=True);b.x=cells
    b.op('live_stain','SourcegreenU87/orange3T3separatestainstocks37C15min; exactdyechemicalattachmentunassembled',S.MAGIC,(b.stock('green_dye','SourceC3100MPphysicalcellstain',K.MOLECULAR,R.STRUCTURAL),b.stock('orange_dye','SourceC34851physicalcellstain',K.MOLECULAR,R.STRUCTURAL)),product='Sourceseparatelystained/dissociatedcellsinTBS; notoneuniversaldual-dyelabeledcell')
    b.clean('microwells','Source1700cells/ulload3minsettle/TBSflushuntrappedcells; retainmicrowellcells')
    b.op('bead_load','Sourceliveimagebefore500beads/ulbeadload/TBSflush',S.MAGIC,((production,R.STRUCTURAL),),product='Sourcecell/beadmicrowells; occupancy/cross-talkunassigned, no perfect1:1assumption')
    b.op('lysis_capture','SourceTCL99percent/mercaptoethanol1percent thenoil seal;50C20min25C90minRNAcapture',S.MAGIC,(b.cond('lysis_stock','SourceTCL/mercaptoethanol/perfluorinatedoilcompartments'),),product='Copresentcaptured/freeRNA/beads/lysateinoil-isolatedmicrowells; exacthistoricalcapturebodyunassembled')
    b.clean('unseal_wash','Sourcewashbuffer/.02UulSUPERaseINunseals/removesuncapturedRNA')
    t=g.product('retained_rna','Co-retainedRNAtemplatecomponent',K.MOLECULAR);primer=g.product('retained_bead_primer','Samephysicallyretainedbead-boundpolyTprimercomponent; fullcommercialsequenceunassembled',K.OLIGO);g.react('capture_components','RepresentcopresentRNAandbeadprimer, nofreshcaptureprimer/intramolecularSELFclaim',S.GENERIC,((b.x,R.TEMPLATE),),(t,primer),changes_molecule=False);b.x=t
    tso=g.input('tso','Source2.5uMTSO; wholeorderingunassembled, noDropseqinheritance',K.OLIGO)
    b.op('rt','SourceMaximaHMinus10Uul/dNTP1mM/SUPERaseIN1Uul/Tween.1percent25C30min42C90min',S.REVERSE_TRANSCRIPTION,((primer,R.STRUCTURAL),(tso,R.STRUCTURAL),b.stock('maxima','SourceMaximaHMinusEP0752'),b.stock('rt_nt','Source1mMdNTPmolecularsubstrate',K.MOLECULAR,R.STRUCTURAL),b.cond('rt_conditions','Source1xMaximaRTbuffer/SUPERaseIN/Tween')))
    protect=[]
    for n,c in enumerate(CODES):protect.append(b.stock('protect'+str(n+1),'SourceunlabeledOBOcomplement1uMeach, distinctfromCy5opticalprobe',K.OLIGO,R.STRUCTURAL,chemical_state={'source_binding_target':c,'whole_unlabeled_ordering':'unassembled; samewholebodyasCy5probe notestablished'}))
    b.op('protect','Source12unlabeledcomplements1uMeach/washbuffer1hRT',S.HYBRIDIZATION,protect)
    b.clean('protect_wash','Sourcewashbufferafterunlabeledcomplements')
    b.op('exoi','SourceExoI1Uul1xbuffer37C45min',S.DIGESTION,(b.stock('exoi_enzyme','SourceExoIM0293L'),))
    b.clean('exoi_wash','SourceTE/TWthen150mMNaOH; retainbeadscDNA/opticalsites')
    for n,c in enumerate(CODES):
        p=g.input('cy5probe'+str(n+1),'Sourceopticalcycle'+str(n+1)+'200nMCy5probe',K.OLIGO,state=oligo_scene(revcomp(c),mod5='Cy5'));b.op('optical'+str(n+1),'Sourcecyclebackgroundscan/probehybridRT10min/wash/probescan',S.HYBRIDIZATION,((p,R.STRUCTURAL),),product='Sourceprobehybridizedbead/cDNAmixture; exactoccupancyunassigned')
        b.clean('optical_strip'+str(n+1),'Source150mMNaOHRT10minstrip/wash; retainbeads/cDNA, fluorescentprobephysicalreleasedfraction')
    region,otherregions=b.split('regions','Sourcecutdeviceinto10separateregions',[('selected','SelectedregionPDMS/beads, indexidentityunassigned',None),('other9','Other9physicallyseparateregions forindependentlibraryprep',None)])
    b.op('region_extract','Source100percentethanol/vortex/sonicate/centrifuge; pelletbeads andremovePDMS',S.MAGIC,(b.cond('extract_conditions','Sourceethanolextractionconditions'),),product='Sourcebeadpellet; residualPDMS/ethanolfractionsseparatedbelow')
    b.clean('region_wash','SourcephysicalPDMSremoval/150mMNaOH/TE-SDS/TE-TW/waterwashes')
    pcrs=b.split('region_two_pcr','Sourceeachregionbeadssplitintotwo50ulPCRs',[('one','First50ulPCRbeadtemplate',None),('two','Second50ulPCRbeadtemplate',None)])
    froms=[]
    for j,templ in enumerate(pcrs):
        b.x=templ;st=g.input('region_pcr_stock'+str(j),'SourcepreparedKAPAHifiHotStart1x/SMRTpcr1uM50ulreaction'+str(j),K.MOLECULAR);cat=g.product('region_cat'+str(j),'SamePCRstockcatalystcomponent',K.CATALYST);nts=g.product('region_nt'+str(j),'SamePCRstocknucleotidecomponent',K.MOLECULAR);pr=g.product('region_primer'+str(j),'SamePCRstockSMRTpcrspecies',K.OLIGO,state=oligo_scene(SMRT));g.react('region_components'+str(j),'RepresentonesourcepreparedPCRreaction, nofreshNTaddition',S.GENERIC,((st,R.TEMPLATE),),(cat,nts,pr),changes_molecule=False)
        b.op('region_pcr'+str(j),'Source12cycles4x98/65/72then8x98/67/72',S.PCR,((cat,R.CATALYST),(nts,R.STRUCTURAL),(pr,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='SourcesingleSMRTspecies, completeinwardsites/productunassembled.',symmetric_species=True));froms.append(b.x)
    rpool=g.product('region_pcr_pool','ActualselectedregiontwoPCRproductpool',K.MOLECULAR);g.react('region_pool','SourcepooltwoPCRproductsofSAMEregion',S.POOL,tuple((x,R.TEMPLATE)for x in froms),(rpool,));b.x=rpool;b.clean('region_spri','SourceSPRI.6:1')
    b.split('tag_aliquot','Source.6ngcDNAselectedpertagmentation',[('selected','Source.6ngselectedcDNA',None),('rest','UntransferredamplifiedregionalcDNAstock',Disposition.UNRESOLVED)])
    b.op('tagment','SourcehistoricalNexteraXTFC1311024.6nginput; loadedstockcompositionunassembled',S.MAGIC,(b.stock('xt_stock','ActualsourceXTtagmentationmolecularpreparation, loadedoligowholeorderingsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.op('regional_index_pcr','Sourcei7regionindex/universalP5REPLACESi5forselectivecDNA5end, wholeindexprimer/cyclesunassembled',S.MAGIC,(b.stock('index_p5','SourceuniversalP5 '+UNIVERSAL_RAW,K.OLIGO,R.STRUCTURAL,state=oligo_scene(UNIVERSAL,mod3='terminal2PSbonds')),b.stock('index_stock','SourceactualNexteraPCRreactionstocki7/enzyme/nt; exactcompositionunassembled',K.MOLECULAR,R.STRUCTURAL)))
    b.clean('index_spri1','SourceSPRI.6:1');b.clean('index_spri2','SourceSPRI1:1');b.split('regional_load','SourcequalifiedregionalRNA-libraryaliquotforpool',[('selected','Selectedregionlibraryaliquot',None),('rest','Untransferredregionlibrarystock',Disposition.UNRESOLVED)])
    otherlib=g.product('other_region_aliquots','Sourceother9separatelypreparedqualifiedRNA-libraryaliquotsforpool',K.MOLECULAR);otherrest=g.product('other_region_remainders','Untransferredotherregionpreparedlibrary/cDNAstocksandphysicalworkupfractions',K.MOLECULAR)
    g.react('other_region_libraries','Other9physicalregionlibraries:beadextraction/2PCR12cycles/0.6ngXT/universalP5+i7/2SPRI; per-regionwholeproductsunassembled',S.MAGIC,((otherregions,R.TEMPLATE),b.stock('other_region_preparation','ActualotherregionseparateRNA-libraryreagentpreparations',K.MOLECULAR,R.STRUCTURAL)),(otherlib,otherrest),reason='Sourceparallelregionpreparationpositive; nocrossregionDNAtemplatepoolinvented.');g.mark_terminal(Terminal(otherrest.id,Disposition.UNRESOLVED,'Otherregionstocks/workupfractions.'))
    pool=g.product('rna_pooled_library','ActualpooledregionalRNA-libraryaliquots',K.MOLECULAR);g.react('rna_pool','SourcepoolseparatequalifiedregionlibrariesAFTERindexing',S.POOL,((b.x,R.TEMPLATE),(otherlib,R.TEMPLATE)),(pool,))
    phix=g.input('rna_phix','Source20percentPhiXDNAcontrollibrary',K.MOLECULAR);out=g.product('rna_library','SourceRNA-derivedlibrarypool+20percentPhiX',K.LIBRARY);g.react('rna_phix_mix','SourceactualPhiXspike',S.GENERIC,((pool,R.TEMPLATE),(phix,R.SPIKE_IN)),(out,),changes_molecule=False);g.sequence(out,SequencingPlan(limitation='SourceNextSeq50021R1/66R2/8I1 customRead1 '+READ1+';12ntcellbarcodeandUMIexistenceexplicit, UMIlength/fullbeadsequence/wholelandingunassembled. TechnicalDNAlookupisnotsecondbiologicalmodality. Optical/celllinking/computationalfiltersnotphysicaldiscard.'))
    g.add_source_item(SourceItem('unassembled','Fullcommercialcaptureoligo/TSO/junctions/readprimerlanding/history-specificstockformulationsunassembled.','Source3originalcapclosed, noabsenceorfamilyrecipeinheritance.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('Sp18isnotnucleotide',OpticalOligo(CODES[0]).nucleotide_count,30);c('AdditionalterminalddCretained',OpticalOligo(CODES[0]).parts[-1],('modified_DNA','C','3ddC'));c('Uretained',OpticalOligo(CODES[0]).prefix.count('U'),1)
    c('LookupDNAreleasedbeforefirstopticaligation',g.steps['release1'].outputs[1],'factory1_supernatants');c('Secondlookupplatephysicallyretained',g.steps['other_lookup127'].inputs[1].material,'factory2_supernatants');c('DNApoolafterseparateindexedPCR',g.steps['lookup_pool128'].inputs[0].material,'lookup_pcr5_selected');c('Retainedlookupprimer',any(p.material=='lookup_retained_primer'and p.role==R.STRUCTURAL for p in g.steps['lookup_extend'].inputs),True)
    c('FactorybeadssupplyRNAcapture',any(p.material=='bead_aliquot_production'for p in g.steps['bead_load'].inputs),True);c('RTusesretainedexternalbeadprimer',any(p.material=='retained_bead_primer'and p.role==R.STRUCTURAL for p in g.steps['rt'].inputs),True);c('NoSELFprimingclaim',any(p.role==R.SELF_PRIMER for p in g.steps['rt'].inputs),False);c('Opticalcycles12',sum('optical_strip'+str(i)in g.steps for i in range(1,13)),12);c('Regionpooltakes2PCRproducts',len(g.steps['region_pool'].inputs),2);c('Otherregionlibrariesderivephysicalpieces',g.steps['other_region_libraries'].inputs[0].material,'regions_other9');c('PostelectrophoresisSYBRstain',g.steps['lookup_stain'].inputs[0].material,'lookup_gel_run_product');c('PhysicalPhiX2branches',sum(any(p.role==R.SPIKE_IN for p in g.steps[k].inputs)for k in ['lookup_phix_mix','rna_phix_mix']),2)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s13059-018-1607-x_supp_13059_2018_1607_MOESM4_ESM.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1186/s13059-018-1607-x')):
        text=re.sub(r'\s+','',p.read_text())
        for n,code in enumerate(CODES):c('OriginalwholemodifiedOBO'+str(n+1),OpticalOligo(code).ordering in text,True);c('OriginalwholeCy5probe'+str(n+1),('/5Cy5/'+revcomp(code))in text,True)
        for name,seq in [('UniversalP5',UNIVERSAL_RAW),('NexterapolyA',NEXTERAPOLYA),('SMRTpcr',SMRT),('customRead1',READ1),('lookupP7primer1',LOOKUP_P7)]:c('Originalwholeprimer'+name,seq in text,True)
    c('Allpartialopticalcodepairscomputed',sum(len(optical_pair(i).strands)for i in range(12)),24);c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
