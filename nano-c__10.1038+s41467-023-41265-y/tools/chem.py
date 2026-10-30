"""Source-qualified private Nano-C partial model; DNA contacts, native RNA readout."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
AUTHORITY='PublishedPMC10497529 main231–245/113–121/Fig3PDFp5; SI1actualFigS4p7; Data4MOESM7 originalXLSX A72/B75–B104. Closed3originals.'
LINKER='TTAGA'
SBFI='CCTGCAGG'
SPACER='TTAGA'
T7='TAATACGACTCACTATAGGGAG'
UP='CTAGGGAGTAAAATACAACAAAGTTGTCAA'
DOWN='ATCGGCTAGAATTTGGCTACACAAACCAGT'
def probe_scene(side='up'):
    if side not in ('up','down'):raise ValueError('Select up/down source probe')
    s=Scene();s.strand('probe',[Segment('First T sourcebiotinylatedbase','T',note='Base biotinylation per Data4A72; terminal5linker chemistry unassigned'),Segment('remaininglinker',LINKER[1:]),Segment('SbfI recognition',SBFI),Segment('spacer',SPACER),Segment('T7 promoter',T7),Segment('30nttargetarm',UP if side=='up'else DOWN)])
    return s
def seq(sc):return ''.join(s.top for st in sc.strands.values()for s in st.segs)
def binding_scene(side='up'):
    s=probe_scene(side);arm=UP if side=='up'else DOWN
    s.anneal('target',[Segment('INFERRED complementary genomic interval',revcomp(arm),inferred=True)],to='probe',pair=('INFERRED complementary genomic interval','30nttargetarm'));return s
def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('Whole molecules unassembled','[unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph('Nano-C Cas9/T7 capture and nativeRNA contact readout',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Actualsource15millioncells; DNAcontactsprimary, exactgenome/allelesunassigned',K.BIOLOGICAL,primary=True)
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,ports=(),kind=S.MAGIC,product='Sourceheterogeneousmixture; wholeproducts/occupancy/offproductsunassembled',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive sourceoperation; fulloccupiedmolecularstates andproductendclassesunassembled, noabsenceclaim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actualretainedfraction; molecularpurity/recoveryunassigned',K.MOLECULAR);w=g.product(k+'_waste','Actualremovedworkupfractions',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Actualphysicalworkup, notperfectspecificity');x=y
    def components(k,label,spec,prepared=False):
        stock=g.input(k,label,K.MOLECULAR);out=[];ports=[]
        for suffix,kind,role,text in spec:
            m=g.product(k+'_'+suffix,text,kind);out.append(m);ports.append((m,role))
        g.react(k+'_components','Represent co-present components of actual prepared stock' if prepared else 'Inventory decomposition only; packaging/numberofstocks/additionschedule unassigned',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ports)
    def allocate(k,label,remainder):
        nonlocal x
        y=g.product(k+'_used',label,K.MOLECULAR);r=g.product(k+'_rest',remainder,K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));x=y
    op('fix','Source2percentformaldehyde10minRT/.125Mglycine',(inp('fix_reagent','Actualformaldehydefixative/glycinequench',K.MOLECULAR,R.STRUCTURAL),))
    clean('lyse','Source2xlysis20minice:10mMTris8/10mMNaCl/.2percentNP40/proteaseinhibitors; actualfractionsunassigned')
    op('sds','Source.5percentSDS62C10min/1percentTritonquench',(inp('sds_supplies','ActualSDS/Tritonchemicalsupplies; packagingunassigned',K.MOLECULAR,R.STRUCTURAL),))
    op('nla_first','SourceNlaIII400U37C4h',(inp('nla_first_stock','Actualfirst400UNlaIII',K.CATALYST,R.CATALYST),))
    op('nla_second','Sourceanother400UNlaIII37Covernight/62C20mininactivation',(inp('nla_second_stock','Actualsecond400UNlaIIIaddition',K.CATALYST,R.CATALYST),))
    op('proximity_ligate','Source5ulHC-T4ligase/1xbuffer/1percentTriton/10ugulBSAasprinted/16C4h',(inp('t4_ligase','Actual5ulHC-T4DNA ligase',K.CATALYST,R.CATALYST),inp('ligation_buffer','Actual1xligationbuffer/Triton'),inp('ligation_bsa','ActualBSA molecularadditiveatprinted10ugul; stockvolume/packagingunassigned, samerecipe notduplicateaddition',K.MOLECULAR,R.STRUCTURAL)))
    op('de_crosslink','Source50ulproteinaseK65Covernight',(inp('proteinase','Actual50ulproteinaseK',K.CATALYST,R.CATALYST),))
    clean('dna_recover','Sourcephenol-chloroform-IAA/ethanol3CDNArecovery')
    allocate('dna_aliquot','Actual20ug3CDNAcuttinginput','Untransferred3CDNAstock; quantity/fateunassigned')
    guides=g.input('guide_pool','ActualpreparedgRNApanel; each3300ng, selectedmembers/scaffoldwholeRNAunassembled; productionoutsidegraph',K.OLIGO)
    cas=g.input('cas9','ActualCas9V3,33pmolperguide',K.MOLECULAR)
    complex_=g.product('rnp','ActualgRNA-Cas9complexmixture; freecomponents/offtargetbindingunassigned',K.CATALYST)
    g.react('rnp_assemble','Source25ulNEBuffer3.1RT10min; Methods<=12gRNA recommendation',S.BINDING,((guides,R.TEMPLATE),(cas,R.STRUCTURAL),inp('rnp_buffer','ActualNEBuffer3.1/25ultotal')),(complex_,))
    op('cas_cut','Source20ugDNA+RNP/waterto250ul37Covernight/65C30min',((complex_,R.CATALYST),inp('cas_conditions','Actualcuttingtemperature/waterallocation; wholesiteendsunassembled')))
    op('rnase','Source1250URNaseIf37C45min/70C20min',(inp('rnase_stock','Actual1250URNaseIf',K.CATALYST,R.CATALYST),))
    clean('cut_clean','Source1xAMPure/elute55ulwater')
    allocate('repair_aliquot','Actual53.5ulcutDNAforrepair','Nominal1.5ulof55uleluatenottransferred; fateunassigned')
    op('repair','Source53.5DNA+6.5FFPEbuffer+2FFPERepairMix20C15min',(inp('repair_mix','Actual2ulNEBNextFFPERepairMix; compositionunassembled',K.MOLECULAR,R.STRUCTURAL),inp('repair_buffer','Actual6.5ulFFPERepairBuffer')))
    clean('repair_clean','Source3xAMPure; elutionvolumeunassigned')
    allocate('extension_aliquot','Actual28.5ulrepairedDNA','Untransferredrepairedstock; total/remainderamountunassigned')
    ps=g.input('probe_mix','Actual10ulbiotinT7probe mix:.25ul40uMeach/water; actualselectedpanelunassigned',K.MOLECULAR)
    pu=g.product('probe_up','Data4B75representativeupstreamcomponent; notassertedwholepanel',K.OLIGO,state=probe_scene('up'),chemical_state={'biotin':'FirstnucleotideT biotinylated perA72, exactlinkageunassigned'})
    pd=g.product('probe_down','Data4B76representativedownstreamcomponent; notassertedwholepanel',K.OLIGO,state=probe_scene('down'),chemical_state={'biotin':'FirstnucleotideT biotinylated perA72, exactlinkageunassigned'})
    po=g.product('other_probes','Actualpanelmembership/otherprobecomponentsunassembled; representativepairisnotassertedcompletepanel',K.MOLECULAR)
    g.react('probe_components','Representknownsourceprobecomponentclassesofoneactual10ulmix; chosenpanelunassigned',S.GENERIC,((ps,R.TEMPLATE),),(pu,pd,po),changes_molecule=False)
    op('probe_extend','Source28.5DNA+10probe-mix+10buffer+1dNTP+1OneTaq=printed50.5ul;95C8min/cool1Cper15sto65/68C5min/4C',((pu,R.STRUCTURAL),(pd,R.STRUCTURAL),(po,R.STRUCTURAL),inp('onetaq','Actual1ulOneTaq',K.CATALYST,R.CATALYST),inp('extension_nts','Actual1ul10mMdNTP',K.MOLECULAR,R.STRUCTURAL),inp('extension_buffer','Actual10ul5xOneTaqbuffer')),product='Sourcepartiallydouble-strandedbiotinT7fused3CDNA plusresidualprobes/unextended/offtargetclasses; wholeproductsunassembled')
    op('dilute','Sourcediluteto200ulwater',kind=S.GENERIC,ports=(inp('dilution_water','Actualwateradditionto200ul'),))
    bound=g.product('captured','ActualC1-boundbiotin-containingDNA/probes; notpuretargetDNA',K.MOLECULAR);flow=g.product('capture_flow','Actualnonretainedcapturefraction',K.MOLECULAR)
    g.react('capture','Source5ulC1beadsperprobe/prewash/resuspend200ulBindingWash; addall/RT30minmix',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),inp('beads','ActualC1streptavidinbeads',K.SELECTION,R.SELECTION),inp('capture_buffer','Actual200ul10mMTris7/1mMEDTA/2MNaCl')),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Actualnonretainedcapturefraction, laterfateunassigned'));x=bound
    clean('bead_wash','Source3x200ul1xBindingWash, then43water+5CutSmart37C20min; retainbead-boundmaterial')
    released=g.product('released_dna','ActualSbfIreleased3CDNA; exactstrandends/promoteroccupancyunassembled',K.MOLECULAR);residual=g.product('residual_beads','Actualresidualbead-boundbiotinfragments/unreleasedclasses; fatesunassigned',K.MOLECULAR)
    g.react('sbfi_release','Source43water+5CutSmart/20USbfI37C20minrelease',S.MAGIC,((x,R.TEMPLATE),inp('sbfi','Actual20USbfI',K.CATALYST,R.CATALYST),inp('sbfi_conditions','Actualsourcewater/CutSmartreleaseconditions')),(released,residual),reason='SourceSbfIcleavesrecognitionsiteandreleasesDNA; fullsealedduplex/endcoordinatesunassembled. Biotin-streptavidinbondnotdeclaredcleaved.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'ResidualbeadmaterialnotreturnedtoIVT'));x=released
    clean('released_clean','Source1xAMPure/elute11ulwater')
    allocate('ivt_aliquot','Actual10ulreleasedDNAforIVT','Nominal1ulof11uleluateuntransferred; fateunassigned')
    ivt_components=components('ivt_mix','Actual3ulT7EnzymeMix, ONEpreparedstock',[
        ('polymerase',K.CATALYST,R.CATALYST,'PositiveT7polymerasecatalystclasscopresentinactual3ulmix; doseunassigned'),
        ('other',K.MOLECULAR,R.STRUCTURAL,'OthercopresentT7EnzymeMixcomponentsunassembled, notfreshadditions')],prepared=True)
    op('ivt','Source10DNA+15twoXT7RiboMAXbuffer+3enzymeMix+2betaine+1inhibitor=31ul/37C60min',ivt_components+(inp('ivt_buffer','Actual15ul2xT7RiboMAXbuffer; stockNTcompositionunassembled',K.MOLECULAR,R.STRUCTURAL),inp('ivt_conditions','Actual2ul5Mbetaine/1ulSUPERaseIn')),product='ActualIVTRNAreportsofDNAcontacts pluscopresentDNA/reagents/offproducts; completeRNAstarts/lengths/yieldunassembled')
    op('poly_a','SourceAM1350PolyATailingKit37C10min',(inp('tailing_supplies','ActualsourcePolyATailingKitreagentinventory; doses/componentsunassembled',K.MOLECULAR,R.STRUCTURAL),),product='SourcepolyAtail-bearingIVTRNApopulation pluscopresenttemplates/reagents/failedproducts; exacttailunassigned')
    clean('rna_clean','Source100ulRNACleanXP/elute12ulwater')
    qc=g.product('rna_qc','ActualQubitRNAHSQCaliquot; volumeunassigned',K.MOLECULAR);prod=g.product('rna_production','RemainingactualRNAproductionstock; volumeunassigned',K.MOLECULAR)
    g.react('rna_qc_split','SeparateactualQubitQCfromproduction',S.SPLIT,((x,R.TEMPLATE),),(qc,prod));g.mark_terminal(Terminal(qc.id,Disposition.QC,'ActualQubitRNAHS'));x=prod
    first_supplies=components('rna002_first','ActualhistoricalRNA002firstligation reagentinventory; packagingunassigned',[
        ('adaptor',K.OLIGO,R.STRUCTURAL,'Sourceknownfirstadaptoroligoclass; wholebody/modifications/strandsunassembled'),
        ('other',K.MOLECULAR,R.STRUCTURAL,'Remainingactualfirstligationreagents; catalyst/NT/bufferidentities/dosesunassembled')])
    op('kit_first_ligation','SourcehistoricalSQK-RNA002firstadaptorligation15min',first_supplies)
    rt_supplies=components('rna002_rt','ActualhistoricalRNA002RTreagentinventory; packagingunassigned',[
        ('enzyme',K.CATALYST,R.CATALYST,'PositiveRTenzymeclass; actualidentity/doseunassembled'),
        ('nts',K.MOLECULAR,R.STRUCTURAL,'PositivecopresentRTnucleotideclass; dose/stockoriginunassembled, notextraaddition'),
        ('other',K.MOLECULAR,R.STRUCTURAL,'RemainingactualRTsupplies; identities/dosesunassembled')])
    op('kit_rt','SourceRNA002RT50C30minplus1ulSUPERaseIn',rt_supplies+(inp('rt_inhibitor','Actualadditional1ulSUPERaseIn'),),product='NativeRNAwithstabilizingcDNA/reagents/uncopiedclasses; RNAremainsporesubstrate, wholeduplexunassembled')
    second_supplies=components('rna002_second','ActualhistoricalRNA002secondligation reagentinventory; packagingunassigned',[
        ('adaptor',K.OLIGO,R.STRUCTURAL,'Sourceknownsecondadaptoroligoclass; wholemotorbearingbody/modifications/strandsunassembled'),
        ('other',K.MOLECULAR,R.STRUCTURAL,'Remainingactualsecondligationreagents; catalyst/NT/bufferidentities/dosesunassembled')])
    op('kit_second_ligation','SourceRNA002secondadaptorligation15min',second_supplies)
    clean('kit_cleanup','Sourcefinal24ulRNACleanXP instead40ul; otherkitcleanup/elutionunassembled')
    allocate('loading','ActualMinIONloadingallocation; volumeunassigned','Untransferredfinallibrarystock; volume/fateunassigned')
    lib=g.product('native_rna_library','ActualNano-CnativeRNAsequencinglibraryreportingDNAcontacts; wholegeometryunassembled',K.LIBRARY)
    g.react('load','SourceFLO-MIN106/R9/MinIONMK2.0/48–72h',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='NativeRNAreadoutofDNAcontacts. HistoricalRNA002wholeadapters/readlanding/tail/featurecyclesunassembled; nocDNAreadoutorRNA004inheritance. Viewpoint/mapping/contactfilterscomputational, notDNAUMIsorphysicalselections.'))
    g.add_source_item(SourceItem('probe_modification','Data4A72 firstbasebiotinylated; B75/B76wholeprobes/T7/SbfI/targetarmsknown.','No5primebiotinlinkerchemistryinvented. binding_scene inferredtargetdefault; occupiedwholeDNAproductsandSbfIendsunassembled.'))
    g.add_source_item(SourceItem('source_quantities','Main231printedBSA10ugul;239recipe50.5ul;243IVT31ul.','Retainactualvalueswithoutnormalizingorpublicationerrorclaim.55→53.5,11→10and28.5fromunknownrepairedstockexplicit.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'s41467-023-41265-y_supp_41467_2023_41265_MOESM7_ESM.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-023-41265-y')):
        t=p.read_text();c('Wholeupstreamsourceprobe',('Upstream probe: '+seq(probe_scene('up')))in t,True);c('Wholedownstreamsourceprobe',('Downstream probe: '+seq(probe_scene('down')))in t,True);c('Actualfirstbasebiotin','The first base in each probe is biotinylated'in t,True)
    for side in ('up','down'):
        sc=binding_scene(side);c('Conditionaltarget '+side,sc.strands['target'].segs[0].inferred,True)
    try:probe_scene('arbitrary');bad=False
    except ValueError:bad=True
    c('Rejectunknownprobeside',bad,True)
    c('SecondNlaIIIisfreshaddition',g.steps['nla_second'].inputs[1].material,'nla_second_stock')
    c('RNaseafterCas9',g.steps['rnase'].inputs[0].material,'cas_cut_product')
    c('Repairuses53.5aliquot',g.steps['repair'].inputs[0].material,'repair_aliquot_used')
    c('Cutstock1.5notwaste',g.terminals['repair_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    c('Extensionuses28.5aliquot',g.steps['probe_extend'].inputs[0].material,'extension_aliquot_used')
    c('Releaseuseswashedbeadfraction',g.steps['sbfi_release'].inputs[0].material,'bead_wash_kept')
    c('ResidualbeadsnotIVT',g.terminals['residual_beads'].disposition,Disposition.UNRESOLVED)
    c('IVTuses10of11',g.steps['ivt'].inputs[0].material,'ivt_aliquot_used')
    c('RNAQCseparate',g.terminals['rna_qc'].disposition,Disposition.QC)
    c('RTbetweenligations',g.steps['kit_rt'].inputs[0].material,'kit_first_ligation_product')
    c('SecondligationafterRT',g.steps['kit_second_ligation'].inputs[0].material,'kit_rt_product')
    c('OnlynativeRNAproductionsequence',set(g.sequencing),{'native_rna_library'})
    c('NoPCRinvented',all(s.kind!=S.PCR for s in g.steps.values()),True)
    c('T7enzymefromsame3ulmix',g.steps['ivt_mix_components'].inputs[0].material,'ivt_mix')
    c('Firstadaptorclassisexplicit',g.materials['rna002_first_adaptor'].kind,K.OLIGO)
    c('Secondadaptorclassisexplicit',g.materials['rna002_second_adaptor'].kind,K.OLIGO)
    c('BSAismolecular',g.materials['ligation_bsa'].kind,K.MOLECULAR)
    c('Wholeproductsunassembled',g.materials['probe_extend_product'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
