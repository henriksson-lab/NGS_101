"""Private source-qualified SIF-seq reporter construction, cell selection and NGS recovery."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,Construct,revcomp
from illumina import STEM
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry,PCRPrimer
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
AUTHORITY='AssociatedAMPMC4008384/EMS57097 Methods113–145/Results82–86; publishedSI19actualFigS1p2/TableS2p9firstsixwholeorders/Notep10–12vectorroles. PublishedmainPDFHTML-notPDF, no fullfinal/AMequivalence. SelectedmNanogBACRP24-73P7/insulatedH19/E14ESundifferentiatedbranch.'
SOURCE_ORDERS={'adaptor_a':'GGGGACAAGTTTGTACAAAAAAGCAGGCTCGCTCTTCCGATC*T','adaptor_b':'/5Phos/GATCGGAAGAGCGACCCAGCTTTCTTGTACAAAGTGGTCCCC','forward':'GGGGACAAGTTTGTACAAAAAAGCAGGCT','reverse':'GGGGACCACTTTGTACAAGAAAGCTGGGT','yfp_f':'GAAGCGCGATCACATGGT','yfp_r':'CCATGCCGAGAGTGATCC'}
CORE_B=STEM+'GA';CORE_A=revcomp(CORE_B)
A=SOURCE_ORDERS['adaptor_a'].replace('*','');B=SOURCE_ORDERS['adaptor_b'].removeprefix('/5Phos/')
F=SOURCE_ORDERS['forward'];REV=SOURCE_ORDERS['reverse']
TRANSFECTIONS=10;RECOVERY_PCRS=5

def oligo_scene(name):
    s=Scene()
    if name=='adaptor_a':s.strand('oligo',[Segment('unpaired_attb1_arm',A[:28]),Segment('paired_core',CORE_A),Segment('terminal_3T','T',note='SourcePSbondimmediatelybeforethisphysicalT; no3primeblock')])
    elif name=='adaptor_b':s.strand('oligo',[Segment('paired_core',CORE_B),Segment('unpaired_attb2_arm',B[14:])],mod5='P')
    else:s.strand('oligo',[Segment(name,SOURCE_ORDERS[name])])
    return s

def adaptor_scene():
    if A!=A[:28]+CORE_A+'T' or B!=CORE_B+B[14:] or len(A)!=43 or len(B)!=42:raise ValueError('SourcewholeYadaptorchanged')
    s=Scene();s.strand('attb1',[Segment('unpaired_attb1_arm',A[:28]),Segment('paired_a',CORE_A),Segment('terminal_3T','T',note='SourcePSbondbeforethisT; ligationendnotblocked')])
    s.anneal('attb2',[Segment('paired_b',CORE_B),Segment('unpaired_attb2_arm',B[14:])],to='attb1',pair=('paired_b','paired_a'),unpaired=('unpaired_attb2_arm',),mod5='P')
    s.note('attb1','Knownwhole43/42ntsourceYadaptor:14bpcoreonly, two28ntarmsunpaired; unpairedtailannotationdoesnotclaimadditionalduplexorfreeoligos.')
    s.ligation_end('attb1','terminal_3T',kind='cohesive',polarity='3-prime',owner='adapter',text='Source3Tprotrusion/attB2opposing5Pcompatiblewithsource3A-tailedBACDNA; wholeoccupied/sealedproductsunassembled')
    return s

def cloning_reference():
    """Conditional two-ended reference; heterogenous genomic molecules are not assigned."""
    return Construct([Segment('left_source_adapter_context',A,inferred=True),Segment('unknown_BAC_insert','[BAC fragment]',placeholder=True,inferred=True,note='Actualfragmentsequence/bothoriginalendsunassigned; intendedligationcontextonly'),Segment('right_source_adapter_context',B,inferred=True)],'INFERRED source-derived intendedcloning-templateclass')

def recovery_reference():
    return Construct([Segment('primer_derived_forward_end',F,inferred=True),Segment('unknown_integrated_amplicon','[integrated insert]',placeholder=True,inferred=True,note='WholeGateway/genomicjunctions notborrowedfromreferencevector'),Segment('primer_derived_reverse_end',revcomp(REV),inferred=True)],'INFERRED primer-derivedampliconends; wholeintegratedinsertunassembled')

def initial_reverse_binding_scene():
    s=Scene();s.strand('adapter_b_context',[Segment('source_b_prefix',B[:13],inferred=True),Segment('reverse_primer_window',B[13:],inferred=True)],mod5='BAC continuation')
    s.anneal('attB2R',[Segment('source_reverse_primer',REV)],to='adapter_b_context',pair=('source_reverse_primer','reverse_primer_window'))
    s.note('adapter_b_context','INFERRED ligatedtemplatecontextonly: wholeknownBsequence/Rprimerinterface, nofreeB5P orknownupstreamBACjunction. InitialYarmsarenotassumedalreadycopiedintofullduplex.')
    return s

def cloning_pcr_geometry():
    return PCRGeometry.unavailable(initial_reverse_binding_scene(),reason='SourcePCRcomplementsYadaptors. KnownlocalattB2Rinterfacecomputed; fullinitialtwo-primer/genomicgeometryunassembled. Do notclaimattB1Fsitealreadyexistsbeforeadaptercompletion.')

def completed_cloning_pcr_geometry():
    return PCRGeometry.exact(cloning_reference(),PCRPrimer('attB1F',F,anneal_nt=len(F)),PCRPrimer('attB2R',REV,anneal_nt=len(REV)),statement='INFERRED intendedcompleted-adaptorPCRreferenceonly: bothsourceprimerlandingscomputedafterYarmscopied. Not the actualinitialPCRsubstrate/first-cyclegeometry.')

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('unknown_integrated_template','[unassembled]',placeholder=True,inferred=True)]);return s

def _same_state(x):
    out={}
    if x.state is not None:out['state']=x.state
    if x.chemical_state is not None:out['chemical_state']=x.chemical_state.copy()
    return out

class Trace:
    def __init__(self,g,x):self.g=g;self.x=x
    def inp(self,key,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.inp(key,label,K.CATALYST,R.CATALYST)
    def op(self,key,label,ports=(),kind=S.MAGIC,product=None,**kw):
        state=kw.pop('state',None);chemical_state=kw.pop('chemical_state',None)
        y=self.g.product(key+'_product',product or 'Actualsourceoperation; wholeheterogeneousmolecules/occupiedjunctions/conversion/offproductsunassembled',K.MOLECULAR,**({'state':state}if state is not None else{}),**({'chemical_state':chemical_state}if chemical_state else{}))
        if kind==S.MAGIC:kw.setdefault('reason','Actualsourceoperation; fulloccupiedsubstrates/products/atommap/yieldunassembled, notabsence.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
    def cleanup(self,key,label,ports=(),preserve=False):
        y=self.g.product(key+'_kept','Actualretainedsourceworkupfraction',K.MOLECULAR,**(_same_state(self.x)if preserve else{}));w=self.g.product(key+'_excluded','Actualexcludedworkupfraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(ports)+(self.inp(key+'_conditions',label),),(y,w));self.g.discard(w,'Sourceexcludedphysicalfraction, recoveryunassigned');self.x=y;return y
    def allocation(self,key,label,selected,rest,qc=False,kind=K.MOLECULAR):
        a=self.g.product(key+'_selected',selected,kind,**_same_state(self.x));r=self.g.product(key+'_remainder',rest,self.x.kind,**_same_state(self.x))
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r));self.g.mark_terminal(Terminal(a.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,selected if qc else rest));self.x=r if qc else a;return self.x
    def components(self,key,label,items):
        stock=self.g.input(key,label,K.MOLECULAR);out=[];ports=[]
        for suffix,kind,role,text in items:
            y=self.g.product(key+'_'+suffix,text,kind);out.append(y);ports.append((y,role))
        self.g.react(key+'_components','Sameactualsuppliesinventorycomponents; stockpackagingunassignedunlessnamedpreparedmix, noextrafreshdose',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ports)

def protocol_graph():
    g=ProtocolGraph('SIF-seq selectedmNanogBAC/insulatedH19/ESYFPselectionandPacBioDNArecovery',authority=AUTHORITY,representation='published_instance')
    bac=g.input('bac_dna','ActualpreparedBACRP24-73P7DNA; entiresequence/shearedfragmentidentitiesunassigned',K.MOLECULAR,primary=True);b=Trace(g,bac)
    g.add_source_item(SourceItem('source_version','PublishedmainFunction-based identification of mammalian enhancers using site-specific integration; actualreadmainisassociatedAMEMS57097/PMC4008384.','FinalPDFreturnedHTML. Wholefinal/AMequivalenceunverified. PublishedSI19sourceoligos/Noterolesusedlocally, no fullvendor/Gateway/SMRTrecipe inherited.'))
    g.add_source_item(SourceItem('known_interface','PublishedTableS2attB1PSbeforeterminalT/attB25P; whole43/42ntYorders,14bpcore,28/28ntunpairedarms/one3T.','Onlyintendedtwo-endedreferencecontextinferredbydefault; exactprimerlandingscomputedlocally. No uniforminsert/ligationpurityorpropagationofadaptorPSontonewPCRbackbones.'))
    g.add_source_item(SourceItem('selection_scope','AMsource10transfections/library and5genomicPCRreactions/sample, actualpoolingaftereacharray.','Templateallocatedbeforeeachrealreaction. HATtargetsintendedHprt-restoredsingle-copyclasses; exactgenomicjunctions/offtargets/copydistributionunassembled. SourceYFP/inputlibrariesnotphysicallypooledtogether.'))
    g.add_source_item(SourceItem('source_literals','AMMethods137prints1percentw/vsaline and0.4umstrainer; actualpublishedmainunread.','Preserveliteralstrainerspecification/executeddeviceunassigned, no40umrepair/publicationerrorclaim. No4percentPFAfromdifferentiatedbranchimportedintoselectedESpreparation.'))
    g.add_source_item(SourceItem('ngs_boundary','AMPCRAmplifyingInsertsandSequencing145positivelystates5PCRproducts/sample pooled/AMPure/PacBioRS.','UpstreamplasmidconstructionandFACSnotsequencing. ActualSMRTfinishing/adapter/readprimer/polymerasedoses/barcodes/readlengthunassembled; do notinventmodernSMRTbellkit.'))
    # Source stocks, with actual PCR-use allocations rather than reuse of one whole bottle.
    primers={}
    names=['clone']+[group+str(i)for group in ('sorted','input')for i in range(1,RECOVERY_PCRS+1)]
    for which in ('forward','reverse'):
        stock=g.input(which+'_stock','Actualsource'+which+'primerstock; totalamount/concentrationunassigned',K.OLIGO,state=oligo_scene(which),chemical_state={'source_ordering':SOURCE_ORDERS[which],'ordered_end_modifications':'unassigned'})
        outs=[]
        for use in names:
            y=g.product(which+'_'+use,'Actualprimerallocationfor'+use+'reaction; doseunassigned',K.OLIGO,**_same_state(stock));outs.append(y);primers[(which,use)]=y
        rest=g.product(which+'_remainder','Untransferredsourceprimerstock; amountunassigned',K.OLIGO,**_same_state(stock));outs.append(rest)
        g.react(which+'_allocations','SourceprimerstockcontributionstoclonePCRand10distinctgenomicPCRreactions',S.SPLIT,((stock,R.TEMPLATE),),tuple(outs));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Exactuntransferredprimeramountunassigned'))
    b.op('shear_bac','SourceBioruptorXLorSonifierII450mechanicalfragmentation',ports=(b.inp('shear_conditions','Actualsourceinstrumentchoiceunassigned'),),kind=S.FRAGMENTATION)
    b.cleanup('size_select','Sourceagarose~1to1.6kb/QIAquickGelExtraction')
    endkit=b.components('end_it_supplies','ActualEnd-Itend-repairkitinventory; exactstockcount/formulationunassigned',[
        ('catalysts',K.CATALYST,R.CATALYST,'Actualsourceend-repaircatalystclass; identitiesunassigned'),('conditions',K.CONDITION,R.CONDITION,'ActualsourceEnd-Itconditions/substratecomponentsunassembled')])
    b.op('end_repair','SourceEnd-ItDNAEnd-Repairkitoperation',ports=endkit,kind=S.END_REPAIR);b.cleanup('repair_cleanup','SourceQIAquickPCRPurification')
    b.op('a_tail','Source50ul/15UKlenowexo-minus/.2mMdATP/1xNEBuffer2/37C30min',ports=(b.cat('klenow','Actual15UKlenow3to5exo-minus'),b.inp('datp','Actual0.2mMdATP',K.MOLECULAR,R.STRUCTURAL),b.inp('tail_conditions','Actual50ul/NEBuffer2/37C30min')),kind=S.TAILING,product='Sourcesource3A-tailedBACfragmentclasses; wholeends/occupancyunassembled');b.cleanup('tail_cleanup','SourceQIAquickPCRPurification')
    oa=g.input('adaptor_a','SourceHPLCAdatpor-attB1withPSbondbeforeterminal3T',K.OLIGO,state=oligo_scene('adaptor_a'),chemical_state={'source_ordering':SOURCE_ORDERS['adaptor_a'],'phosphorothioate_bond_between_zero_based':[41,42],'terminal_nucleotide':'physicalT; not3block'})
    ob=g.input('adaptor_b','SourceHPLCAdaptor-attB2/5P',K.OLIGO,state=oligo_scene('adaptor_b'),chemical_state={'source_ordering':SOURCE_ORDERS['adaptor_b'],'five_prime':'P'})
    ad=g.product('annealed_adapter','Actualequimolar95C/slowcoolYadaptorstock',K.OLIGO,state=adaptor_scene(),chemical_state={'source_whole_orders':{k:SOURCE_ORDERS[k]for k in ('adaptor_a','adaptor_b')},'core_bp':14,'adaptor_a_PS_bond':[41,42],'adaptor_b_5prime':'P'})
    g.react('anneal_adapter','SourceequimolarHPLCstrands95Cheatblock/slowreturnRT',S.BINDING,((oa,R.TEMPLATE),(ob,R.STRUCTURAL),b.inp('anneal_conditions','Actualsourceannealingprogramme')),(ad,))
    tmp=b.x;b.x=ad;ad_use=b.allocation('adapter_use','Actual2uMadaptorcontributionto50ulligation; stockconcentration/volumeunassigned','Actualsourceadaptorligationallocation','ActualuntransferredYadaptorstock',kind=K.OLIGO);b.x=tmp
    lig=b.components('ligation_supplies','ActualNEBNextQuickLigationModule supplies aggregateinventory; notprovenonepremix',[
        ('ligase',K.CATALYST,R.CATALYST,'ActualsourceDNAligaseclass; preciseidentity/stockdoseunassigned'),('buffer',K.CONDITION,R.CONDITION,'ActualsourceQuickLigationconditions/formulationunassembled')])
    b.op('adapter_ligate','Source2uMYadaptor/NEBNextQuickLigationModule',ports=((ad_use,R.STRUCTURAL),)+lig,kind=S.LIGATION,geometry=AdapterLigationGeometry.known('cohesive',adaptor_scene(),statement='KnownsourceY3T/opposing5P andsourceBAC3Aclass; fulljunction/yield/occupiedgenomicendsunassembled'),chemical_state={'source_whole_adaptor_orders':{k:SOURCE_ORDERS[k]for k in ('adaptor_a','adaptor_b')},'reference_scope':'ActualY-endedfragmentclasses, initialarmsnotassumedalreadyfilled; wholeinsert/offproductsunassembled'})
    b.cleanup('ligation_cleanup','SelectedmNanogsource35ulAMPureXPper50ulligation',preserve=True)
    b.allocation('qubit','ActualQubitdsDNAHSanalyticalaliquot; amountunassigned','ActualQubitQCcontribution','ActualpostQCDNAstock',qc=True)
    b.allocation('clone_pcr_input','Actual5to10nglibraryDNAfor50ulPCR','Actual5to10ngPCRtemplateclass','UntransferredadaptedDNAstock')
    pfu=b.components('pfu_supplies','Actual50ulPfuUltraIIFusionHSPCRsuppliesaggregate; stockpackagingunassigned',[
        ('enzyme',K.CATALYST,R.CATALYST,'ActualPfuUltraIIFusionHSenzymeclass'),('nts',K.MOLECULAR,R.STRUCTURAL,'ActualPCRnucleotideclass; sourceindividualdosesunassigned'),('buffer',K.CONDITION,R.CONDITION,'ActualPCRbufferco-components')])
    b.op('clone_pcr','Source95C2min/20x95C20s65C20s72C15s/72C3min',ports=((primers[('forward','clone')],R.STRUCTURAL),(primers[('reverse','clone')],R.STRUCTURAL))+pfu,kind=S.PCR,geometry=cloning_pcr_geometry(),state=cloning_reference(),chemical_state={'reference_scope':'INFERRED intendedPCRclass; wholeinsert/offproductsunassigned','PCR_copy_modifications':'No automaticinheritanceoforiginaladaptorPSor5Pontoallnewstrands'})
    b.cleanup('clone_size','Source1percentagarose1to1.6kb/QIAquickelute25ul',preserve=True)
    b.allocation('bp_input','Actualupto100ngof25uleluateforBP','Actualupto100ngampliconlibrary','Untransferred25ulDNAeluatestock')
    donor=g.input('pdonr221','Actual200ngpDONR221; wholevectorunassembled',K.MOLECULAR)
    b.op('gateway_bp','Source15ulBP/3ulBPClonaseII/200ngpDONR221/RT~20h',ports=((donor,R.TEMPLATE),b.cat('bp_clonase','Actual3ulBPClonaseIIpreparedenzymemix'),b.inp('bp_conditions','Actual15ulreactionRT~20h')),kind=S.SITE_SPECIFIC_RECOMBINATION,product='SourceBPmixturewithentryproduct/unreactedtemplates/byproducts; exactjunctionsunassembled')
    b.allocation('bp_transfer','Actual10of15ulBPmixturetoLR','Actual10ulBPforLR','Actualnominal5ulBPstockuntransferred')
    dest=g.input('h19_vector','Actual2ul150ngulfrozenpreparedpSKB1-Venus-H19destinationvector; exactfullmapunassembled',K.MOLECULAR,chemical_state={'source_roles':['Hprthomologyarms','Gatewaycloningsite','minimalhsp68promoter','VenusYFP','H19insulatorbeforeHprtpromoter'],'full_sequence_scope':'PresentinselectedSIbutnotfullyassembled/coordinatevalidatedinthispass'})
    b.op('gateway_lr','Source10ulBP+2ulvector+3ulLRClonaseII/RT~20h',ports=((dest,R.TEMPLATE),b.cat('lr_clonase','Actual3ulLRClonaseIIpreparedenzymemix'),b.inp('lr_conditions','SourceRT~20h; listednominal15ultotal')),kind=S.SITE_SPECIFIC_RECOMBINATION,product='SourceLRmixturewithreporterplasmidlibrary/unreactedDNA/byproducts; completevectorsunassembled')
    b.allocation('lr_transform','Actual1ulofLRmixturefortransformation','Actual1ulLRDNAcontribution','UntransferredLRmixture; nominal14uliflistedvolumesadd')
    bacteria=g.input('top10_cells','Actual50ulOneShotTOP10electrocompetentEcoli',K.BIOLOGICAL)
    b.op('transform','Sourceelectroporate1ulLRinto50ulTOP10; SOCrecovery30min',ports=((bacteria,R.TEMPLATE),b.inp('soc','ActualrichSOCmedium contribution, amountunassigned',K.MOLECULAR,R.STRUCTURAL),b.inp('transform_conditions','Actualsourceelectroporationprogrammeunassigned')))
    b.op('bacterial_growth','Source500mlLB/200ugmlampicillin37C~12.5h',ports=(b.inp('lb','Actual500mlLBgrowthmedium',K.MOLECULAR,R.STRUCTURAL),b.inp('ampicillin','Actual200ugmlampicillinselectiondrug',K.MOLECULAR,R.STRUCTURAL),b.inp('growth_conditions','Actual37C~12.5h')),product='Sourcegrownplasmidlibrary/bacterialpopulation; copyyield/clonecompositionunassembled')
    b.cleanup('plasmid_maxi','SourceQiagenPlasmidMaxi/large-low-copyprotocol; numberofpreps/yieldunassigned',ports=(b.inp('maxi_supplies','ActualMaxikitworkupinventoryunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.op('pmei','SourcePmeIplasmidlinearization',ports=(b.cat('pmei_enzyme','ActualsourcePmeI'),b.inp('pmei_conditions','Exactsourcewholevectorcutpositions/doseunassembled')),kind=S.DIGESTION,product='SourcePmeI-linearizedreporterplasmidlibrarycollection; wholeDNAends/junctionsunassembled')
    b.cleanup('linear_dna','Sourcephenol/chloroform/70percentEtOH/recoverTris.1mMEDTA7.5')
    # Ten independently allocated templates/cell contributions, then actual pooled populations.
    templates=[]
    for i in range(1,TRANSFECTIONS+1):templates.append(g.product('transfection_dna_'+str(i),'Actual~20uglinearplasmidDNAallocationforoneof10sourceelectroporations',K.MOLECULAR))
    remainder=g.product('linear_dna_remainder','Untransferredlinearizedplasmidlibrarycollection; totalyield/amountunassigned',K.MOLECULAR)
    g.react('transfection_dna_allocations','AllocateDNA BEFOREtenactual~20ugelectroporations; noyieldassumption',S.SPLIT,((b.x,R.TEMPLATE),),tuple(templates)+(remainder,));g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'SourceuntransferredlinearDNAstockquantityunassigned'))
    selected=[]
    for i,dna in enumerate(templates,1):
        tag='es'+str(i);host=g.input(tag+'_cells','Actual1to1.5x10^7E14Tg2a.4malemouseEScellcontribution; Hprt36kbdeletionexons1–2',K.BIOLOGICAL)
        context=g.product(tag+'_context','Sameactualhostcell/genomiccontextcomponent',K.BIOLOGICAL);hr=g.product(tag+'_endogenous_hr','Endogenoussourcehosthomologousrecombinationmachineryclass; notfreshpurifiedenzyme',K.CATALYST)
        g.react(tag+'_components','Samehostcellinventoryviews, nophysicalisolationoraddedHRenzyme',S.GENERIC,((host,R.TEMPLATE),),(context,hr),changes_molecule=False)
        cells=g.product(tag+'_transfected','Sourceelectroporated/intendedHprt-integratedreporterpopulation; exactjunction/copy/offtargetsunassembled',K.BIOLOGICAL)
        g.react(tag+'_integration','Source~20ugDNA/.8mlHBS/250V500uF/transfectionandhosthomologousintegration',S.INTEGRATION,((dna,R.TEMPLATE),(context,R.TEMPLATE),(hr,R.CATALYST),b.inp(tag+'_hbs','Actual0.8mlHEPESbufferedsaline'),b.inp(tag+'_electroporation','Actual250V500uF')),(cells,))
        hat=g.input(tag+'_hat','ActualsourceHATsupplementedESselectionmedium contribution; doseunassigned',K.SELECTION,chemical_state={'named_source_components':['hypoxanthine','aminopterin','thymidine'],'source_scope':'NamedHATsupplement, no extraseparatefreshnucleosidedosesinvented'})
        yes=g.product(tag+'_selected','ActualHAT-selectedcellpopulation; intendedHprt-restorationenrichmentnotperfectgenotyping',K.BIOLOGICAL);no=g.product(tag+'_excluded','Actualdrug-selectionexcludedcell/materialfraction',K.WASTE)
        g.react(tag+'_hat_select','SourceHATselectionstarts24haftertransfection/3to10days',S.SELECTION,((cells,R.TEMPLATE),(hat,R.SELECTION),b.inp(tag+'_hat_conditions','Actual24hstart/3to10days')),(yes,no));g.discard(no,'Actualsourceexcludeddrug-selectionfraction')
        ht=g.input(tag+'_ht','Actual1xHTsupplementedESmediumcontributionfor2days; doseunassigned',K.MOLECULAR,chemical_state={'source_components':['hypoxanthine','thymidine'],'stockpackaging':'ActualnamedHTsupplement, individualdosesunassigned'})
        fed=g.product(tag+'_ht_cells','ActualpostselectionHT-fedsourcecellpopulation',K.BIOLOGICAL)
        g.react(tag+'_ht_feed','SourceHTmedium2daysafterHATselection',S.MAGIC,((yes,R.TEMPLATE),(ht,R.STRUCTURAL),b.inp(tag+'_ht_conditions','Actualsource2days')),(fed,),reason='Actualcellfeeding/growthboundary; completeDNA/cell/mediumchemistry/yieldunassembled.');selected.append(fed)
    pooled=g.product('pooled_es_cells','Actualsame-librarypooloftenindependentselectedtransfectionpopulations',K.BIOLOGICAL)
    g.react('pool_es','Sourcepoolsame-libraryselectedcells BEFOREexpansion/sorting',S.POOL,tuple((x,R.TEMPLATE)for x in selected),(pooled,));b.x=pooled
    b.op('expand_cells','Sourceexpandpooledcells onfreshplates',ports=(b.inp('expansion_medium','ActualfreshESmedium contribution, amountunassigned',K.MOLECULAR,R.STRUCTURAL),b.inp('expansion_conditions','Sourceexpansionprogrammeunassigned')),product='ActualexpandedESreporter-librarypopulation; wholecellnumbers/compositionunassembled')
    unsorted=g.product('input_cells','Actualpre-sortunsortedcellaliquotforinputcontrol',K.BIOLOGICAL);sortcells=g.product('sort_cells','ActualsourcecellcontributionforphysicalYFPsorting',K.BIOLOGICAL)
    g.react('input_control_allocation','SI1actualaliquotcollectedBEFOREsorting, remainingcontributionforsorting; amountsunassigned',S.SPLIT,((b.x,R.TEMPLATE),),(unsorted,sortcells));b.x=sortcells
    b.cleanup('trypsin','SourcePBSwash/trypsinharvest/pellet/removetrypsin/PBSwash',ports=(b.inp('trypsin_reagent','Actualtrypsinharvestcontribution',K.CATALYST,R.CATALYST),))
    b.op('sort_suspension','Source1percentwvsaline/pipette/literal0.4umstrainer',ports=(b.inp('saline','Actualsource1percentwvsaline'),b.inp('strainer','Sourceliteral0.4umstrainer/executedspecificationunassigned')),kind=S.GENERIC)
    yfp=g.product('yfp_cells','ActualFACSselectedrobustYFP-expressingcellpopulation',K.BIOLOGICAL);neg=g.product('nonselected_cells','ActualnonselectedFACSoutcome; amount/fateunassigned',K.BIOLOGICAL)
    g.react('facs','SourceInflux/Spigot6.1.10physicalYFPsorting; controlcalibrationoutsideproductiontrace',S.PARTITION,((b.x,R.TEMPLATE),b.inp('facs_conditions','Actualsourcegatecalibratedpositive/negativecontrols, numericgateunassigned')),(yfp,neg));g.mark_terminal(Terminal(neg.id,Disposition.UNRESOLVED,'NonselectedcellfractionnotNGSinputforselectedsortedbranch'))
    for group,population in [('sorted',yfp),('input',unsorted)]:
        b.x=population;b.cleanup(group+'_gdna','SourceQIAampDNAMini genomicDNAisolation',ports=(b.inp(group+'_gdna_supplies','ActualQIAampDNAworkupinventory; recipevolumesunassembled',K.MOLECULAR,R.STRUCTURAL),))
        aliquots=[g.product(group+'_dna_'+str(i),'Actualupto100nggDNAallocationforoneof5source50ulPCRreactions',K.MOLECULAR)for i in range(1,RECOVERY_PCRS+1)]
        rest=g.product(group+'_gdna_remainder','UntransferredsourcegDNAstock; quantityunassigned',K.MOLECULAR)
        g.react(group+'_dna_allocations','Allocate actualgDNA BEFOREfiveindependentPCRreactions',S.SPLIT,((b.x,R.TEMPLATE),),tuple(aliquots)+(rest,));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'ActualgDNAstockremainderunassigned'));products=[]
        for i,dna in enumerate(aliquots,1):
            tag=group+str(i);sup=b.components(tag+'_pcr_supplies','Actual50ulPlatinumTaqHFPCRsuppliesaggregate; packaging/stockcountunassigned',[
                ('enzyme',K.CATALYST,R.CATALYST,'ActualPlatinumTaqHighFidelityenzymeclass'),('nts',K.MOLECULAR,R.STRUCTURAL,'ActualPCRnucleotidecomponentclass; dosesunassigned'),('buffer',K.CONDITION,R.CONDITION,'ActualPCRbufferco-components')])
            product=g.product(tag+'_pcr_product','ActualindependentrecoveryPCRproductclass; primer-endviewinferred, wholegenomicinsert/offproductsunassembled',K.MOLECULAR,state=recovery_reference(),chemical_state={'reference_scope':'INFERRED primer-derivedDNAends, notwholeintegratedtemplate','sample_group':group})
            g.react(tag+'_pcr','Source50ul/up100nggDNA/94C2min/30x94C30s55C30s68C90s',S.PCR,((dna,R.TEMPLATE),(primers[('forward',tag)],R.STRUCTURAL),(primers[('reverse',tag)],R.STRUCTURAL))+sup,(product,),geometry=PCRGeometry.unavailable(unknown_scene(),reason='KnownwholeattBprimers/sourceactualgenomicPCR; integratedGateway/genomictemplatejunctions/fulloccupiedsitesnotassembled, no vectorrecipeinheritance.'));products.append(product)
        pool=g.product(group+'_pcr_pool','Actualphysicalpoolof5independentsame-samplePCRproducts',K.MOLECULAR,state=recovery_reference(),chemical_state={'reference_scope':'INFERRED primer-derivedends ofreportedampliconclass; fullmixtureunassembled','sample_group':group})
        g.react(group+'_pool_pcr','SourcepoolfiveactualPCRproductsBEFOREAMPurepurification',S.POOL,tuple((x,R.TEMPLATE)for x in products),(pool,));b.x=pool;b.cleanup(group+'_ampure','SourceAMPureXPpurification; ratio/elutionvolumeunassigned',preserve=True)
        b.allocation(group+'_loading','ActualPacBioRSfinishing/loadingcontribution; amountunassigned','Actualsourceampliconallocationforsequencing','Untransferredsampleampliconstock')
        finishing=g.input(group+'_ngs_supplies','ActualreportedPacBioRSfinishing/loadinginventory; precisechemistry/SMRTbell/primer/polymerasedosesunassembled',K.MOLECULAR)
        library=g.product(group+'_library','ActualsamplePacBioRSsequencinglibraryclass; fullmolecules/readgeometryunassembled',K.LIBRARY)
        g.react(group+'_ngs_finish','AMsourcePCRampliconssequencedonPacBioRS; finishingunassembled',S.MAGIC,((b.x,R.TEMPLATE),(finishing,R.STRUCTURAL)),(library,),reason='ActualNGSendpointreported; no modernSMRTbell/adapter/primer/lengthdefaultimported. Upstreamconstruct/FACSstagesnotsequencing.')
        g.sequence(library,SequencingPlan(limitation='ActualsourcePacBioRSampliconsequencing; exactSMRTfinishing/runlength/primerlanding/barcodesunassembled. Sortedversusinputcomparisoncomputational, no physicalmixingofthese2libraries.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('SourceadaptorAphysical43nt',len(A),43);c('SourceadaptorBphysical42nt',len(B),42)
    c('Exact14bpcore',CORE_A,revcomp(CORE_B));adaptor_scene();cloning_pcr_geometry();completed_cloning_pcr_geometry()
    c('Sourceone3Toverhang',A[-1],'T');c('PSbeforeterminalT',SOURCE_ORDERS['adaptor_a'][-3:],'C*T')
    c('Knownprimerbodies29each',(len(F),len(REV)),(29,29))
    c('SourceinitialYdoesnotinventtwo-primercompletedduplex',cloning_pcr_geometry().availability,'unavailable')
    c('Cloningreferenceuncertainbydefault',all(s.inferred for s in cloning_reference()),True)
    c('Recoveryreferenceuncertainbydefault',all(s.inferred for s in recovery_reference()),True)
    c('Actualadapterallocationpreservesstate',g.materials['adapter_use_selected'].state is not None,True)
    c('Ligationusesactualadapterallocation',g.steps['adapter_ligate'].inputs[1].material,'adapter_use_selected')
    c('ClonePCRusesselected5to10ng',g.steps['clone_pcr'].inputs[0].material,'clone_pcr_input_selected')
    c('LRusesactual10of15BP',g.steps['gateway_lr'].inputs[0].material,'bp_transfer_selected')
    c('Transformationusesactual1ulLR',g.steps['transform'].inputs[0].material,'lr_transform_selected')
    c('TenindependenttemplatesallocatedBEFOREtransfection',len(g.steps['transfection_dna_allocations'].outputs)-1,10)
    c('Tentransfectionproductsphysicallypooled',len(g.steps['pool_es'].inputs),10)
    c('FirstintegrationusesactualDNAallocation',g.steps['es1_integration'].inputs[0].material,'transfection_dna_1')
    c('LastintegrationusesdistinctactualDNAallocation',g.steps['es10_integration'].inputs[0].material,'transfection_dna_10')
    c('UnsortedcontrolallocatedbeforeFACS',g.steps['input_control_allocation'].outputs,('input_cells','sort_cells'))
    for group in ('sorted','input'):
        c(group+'fiveDNAtemplatesallocatedbeforePCR',len(g.steps[group+'_dna_allocations'].outputs)-1,5)
        c(group+'actual5PCRproductsinpool',tuple(x.material for x in g.steps[group+'_pool_pcr'].inputs),tuple(group+str(i)+'_pcr_product'for i in range(1,6)))
        c(group+'pooledPCRstockpreservesprimerendview',g.materials[group+'_ampure_kept'].state is not None,True)
        c(group+'NGSusesactualallocation',g.steps[group+'_ngs_finish'].inputs[0].material,group+'_loading_selected')
    c('Nofakeinputsortedlibrarypool',any({x.material for x in s.inputs}=={'sorted_library','input_library'}for s in g.steps.values()),False)
    p=sources_root()/HERE.parent.name/'nmeth.2886_41592_2014_BFnmeth2886_MOESM19_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/nmeth.2886')):
        t=p.read_text()
        for key,raw in SOURCE_ORDERS.items():c('WholepublishedTable2ordering '+key,raw in t,True)
    main=sources_root()/HERE.parent.name/'nmeth.2886_PMC4008384.xml.txt'
    if have(c,Source(main,'https://doi.org/10.1038/nmeth.2886')):
        t=main.read_text();c('SourcefivePCRandpool','Five PCR reactions were performed for each sample and pooled'in t,True);c('Sourcetentransfections','Ten such transfections were performed for each library'in t,True);c('SourceactualPacBioRSendpoint','PCR amplicons were sequenced using a PacBio RS'in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
