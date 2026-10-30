"""Private source-owned endogenous-DPC preparation; published main319 instance."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
TITLE='DPC-seq endogenous crosslinks, technical-preparation component'
AUTHORITY='PublishedOA PMC11098753 main181–187/319;publishedTable6MOESM8originalwholeapplications/orderingsread. Twooriginals; citedProtocolExchange10.21203/rs.3.pex-2562/v1actuallyblocked, no version substitution.'
def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('unassembled','[historicalwholemoleculegeometryunassembled]',placeholder=True)]);return s
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','SourceRPE1technicalreplicate,seeded160000cells/welland24h1percentFBSsynchronization; actualharvestcountunassigned',K.BIOLOGICAL,primary=True)
    def stock(key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(key,label,kind,**kw),role)
    def cond(key,label):return stock(key,label,K.CONDITION,R.CONDITION)
    def op(key,label,kind,inputs=(),product='Sourcephysicalmixturewithresidualsubstrates/reagents; exactgeometryandcompletionunassigned',**kw):
        nonlocal x
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Sourcehistoricalkitcomponent/internalworkupandwholeproductgeometryunassembled.'
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);a=g.product(key+'_kept','Physicallyretainedsourcefraction;recovery/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Sourcephysicallyremovedworkupfraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(a,w));g.discard(w,'Sourcephysicalworkup,noperfectmolecularselectivity.');x=a
    def aliquot(key,label,rest):
        nonlocal x
        a=g.product(key+'_used',label,K.MOLECULAR);b=g.product(key+'_rest',rest,K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,Disposition.UNRESOLVED,rest));x=a
    def qc(key,label):
        nonlocal x
        a=g.product(key+'_kept','ProductionpreparationafterseparateQCaliquot',K.MOLECULAR);q=g.product(key+'_sample',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,q));g.mark_terminal(Terminal(q.id,Disposition.QC,label));x=a
    op('induce','Source1.75mMformaldehyde1h;selectedimmediatecollectionnot6hrecovery',S.CROSS_LINKING,(cond('formaldehyde','Sourceformaldehyde'),),product='SourceFA-treatedcellswithcovalentDPCsandotherlesions;nohomogeneouscrosslinkformationassumed')
    clean('pbs_wash','SourcetwicePBSwashandcellcollection')
    op('lysis','Source150ul2percentSDS/20mMTris7.5scrapelysis;freezeLN2/store-80C/thaw55C5min1200rpm',S.PERMEABILIZATION,(cond('lysis_condition','SourceSDS/Trisandtemperatureprogramme'),),product='Sourcemixeddenaturinglysate;solubleandDPC/proteinmaterialnotyetphysicallypartitioned')
    op('sonicate','SourceCovarisE220evo130ultube,one120scycle;nominalvesselcapacitynotdeclaredsamplealiquot',S.FRAGMENTATION,(cond('sonication_condition','SourceCovarisprogramme'),))
    op('dilute','Source1.5mltubetransfer+270ul2percentSDS/20mMTris7.5',S.GENERIC,(cond('dilution','SourceSDS/Trisaddition'),))
    control=g.product('input_10percent','Source10percenttotallysatephysicalInputaliquot',K.MOLECULAR);production=g.product('production_90percent','Sourceproductionlysateremainder',K.MOLECULAR);g.react('input_split','Physical10percentInputvsproductionremainder',S.SPLIT,((x,R.TEMPLATE),),(control,production));x=production
    ic=g.input('genejet_condition','SourceGeneJETGenomicDNAPurificationKit',K.CONDITION);inputdna=g.product('input_dna','SourceextractedInputDNAcontrol;laterlibrary/preparationpathunassembled',K.MOLECULAR);iw=g.product('input_waste','SourceGeneJETworkupfraction',K.WASTE);g.react('input_extract','Actual10percentInputDNAextraction,notallQCorsequencedbydefault',S.CLEANUP,((control,R.TEMPLATE),(ic,R.CONDITION)),(inputdna,iw));g.discard(iw,'Sourcekitphysicalworkup.');g.mark_terminal(Terminal(inputdna.id,Disposition.UNRESOLVED,'ExtractedInputcontrolisrealphysicalmaterial; laterlibraryscopeunassigned.'))
    pc=g.input('precipitate_condition','Source400ul200mMKCl/20mMTris7.5,ice5min/fullspeed4C5min',K.CONDITION);pellet=g.product('dpc_pellet','Sourceprotein/DPC-richpelletwithpossibleunlinkedDNAcarryover',K.MOLECULAR);soluble=g.product('soluble_dna','Sourcesoluble-DNAsupernatanttransferredforquantification',K.MOLECULAR);g.react('precipitate','SourceKCl-SDSpartitionretainsBOTHobservedfractions',S.CLEANUP,((x,R.TEMPLATE),(pc,R.CONDITION)),(pellet,soluble));x=pellet
    sr=g.input('soluble_rnase_enzyme','Source0.2mgmlDNase-freeRNaseAinSOLUBLEtube',K.CATALYST);sc=g.input('soluble_rnase_condition','Source37C30min',K.CONDITION);sd=g.product('soluble_rnase_product','SourceRNase-treatedsolubleDNApreparationwithRNAfragments',K.MOLECULAR);g.react('soluble_rnase','SourceindependentsolublefractionRNasetreatment',S.DIGESTION,((soluble,R.TEMPLATE),(sr,R.CATALYST),(sc,R.CONDITION)),(sd,))
    sq=g.product('soluble_qubit','SourceQubitHSanalyticalaliquotfromsolubleDNA;allocationunassigned',K.QC);srest=g.product('soluble_remainder','Sourceunmeasuredsolublepreparationremainder',K.MOLECULAR);g.react('soluble_qc','SeparateactualQubitaliquotfromentiresolublestock',S.SPLIT,((sd,R.TEMPLATE),),(sq,srest));g.mark_terminal(Terminal(sq.id,Disposition.QC,'SourcequantificationforDPCratio'));g.mark_terminal(Terminal(srest.id,Disposition.UNRESOLVED,'SourceunmeasuredsolubleDNAisnotimplicitlydiscarded/returnedtoproduction.'))
    for i in range(1,4):clean('wash'+str(i),'Sourcewash'+str(i)+':400ulKClbuffer55C5min1200rpm→ice5min→fullspeed4C5min;retainreprecipitatedpellet')
    op('proteinase','Source400ulKClbuffer+ProteinaseK0.2mgml55C45min800rpm',S.DEPROTEINIZATION,(stock('proteinase_k','SourceProteinaseK'),),product='SourceproteolyzedDPC-richpreparationwithDNA/proteinfragments/residualconjugates;notprovedchemicalreversalofthecrosslink')
    op('bsa','Sourceadd10ulUltraPureBSA,ice5min',S.GENERIC,(stock('bsa_stock','SourceBSAphysicalproteincontribution',K.MOLECULAR,R.STRUCTURAL),))
    cold=g.input('cold_spin_condition','Sourcefullspeed4C5min',K.CONDITION);dna=g.product('released_supernatant','SourcecollectedDNA-containingsupernatantafterproteinprecipitation',K.MOLECULAR);res=g.product('post_pk_pellet','Sourceco-separatedresidualprecipitate',K.MOLECULAR);g.react('post_pk_spin','SourcephysicalsupernatantcollectionafterBSA/coldspin',S.CLEANUP,((x,R.TEMPLATE),(cold,R.CONDITION)),(dna,res));g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Sourcecollectssupernatant;remainingpelletlaterfateunassigned.'));x=dna
    op('rnase','Source0.2mgmlDNase-freeRNaseA37C30mininCROSSLINKEDfraction',S.DIGESTION,(stock('rnase_a','SourcefreshRNaseAinsourceDPCtube'),))
    qc('qubit','SourceQubitHSaliquotfromcrosslinkedDNA;ratioanalysisnotphysicalfractionation')
    aliquot('library_dna','Source50ngDPC-fractionDNAselectedforlibraryconstruction','UntransferredcrosslinkedDNAstock;amount/laterfateunassigned')
    op('precipitate_dna1','Source300mMsodiumacetate+1ulglycogen+2.5xethanol',S.GENERIC,(stock('glycogen1','Source1ulglycogenmolecularcarrier',K.MOLECULAR,R.STRUCTURAL),cond('ethanol1','SourceNaAc/ethanolprecipitationconditions')))
    clean('ethanol1_clean','Sourceethanolprecipitationphysicalrecovery/resuspend20ulnuclease-freewater;carrier/carryoverunassigned')
    op('gel_stain','Source1percentagarosegelstainedSYBR-Gold',S.GENERIC,(stock('sybr_gold','SourceSYBR-GoldDNAbindingdyemolecularcontribution;amountunassigned',K.MOLECULAR,R.STRUCTURAL),cond('gel_condition','Sourceagaroseelectrophoresis/staining')),product='Sourcegel-separated/stainedDNAwithbound/freedyestatesunassembled;noassumedquantitativedyeremoval')
    aliquot('gel_band','Sourceexcised400–1000bpDNA-smeargelslices','Unselectedgelregions/material;laterfateunassigned')
    op('gel_elute','Source500ul10mMTris8/1mMEDTA/0.02percentSDS,rotate4Covernight',S.GENERIC,(cond('gel_extraction_buffer','Sourcegel-extractionbuffer'),),product='Sourcegelslices/bufferwithdiffusedandresidualDNA/dyecarrierclasses')
    fc=g.input('spinx_condition','SourceSpin-X14000g10min4C',K.CONDITION);flow=g.product('spinx_eluate','SourcecollectedDNA-containingcolumneluate',K.MOLECULAR);gel=g.product('spinx_residual','Sourcecolumn/gelresidualmaterial',K.MOLECULAR);g.react('spinx','PhysicalSpin-Xeluateseparation',S.CLEANUP,((x,R.TEMPLATE),(fc,R.CONDITION)),(flow,gel));g.mark_terminal(Terminal(gel.id,Disposition.UNRESOLVED,'Residualgel/columnfractionnotassumedfullyDNA-free.'));x=flow
    op('precipitate_dna2','Sourceethanolprecipitateasbefore;300mMNaAc/1ulglycogen/2.5xethanol',S.GENERIC,(stock('glycogen2','Sourcefreshsecond1ulglycogencarrier',K.MOLECULAR,R.STRUCTURAL),cond('ethanol2','SourceNaAc/ethanolconditions')))
    clean('ethanol2_clean','SourcephysicalDNAprecipitation/resuspend50ulwater')
    op('kit_pre_pcr','SourceNEBNextUltraIIDNAE7645Lpre-PCRpreparation;adapterdilution1/10;wholehistoricalkitstages/endsunassembled',S.MAGIC,(stock('kit_adapter','SourceactuallyselectedhistoricalNEBNextadapteraliquotat1/10stockdilution;ordering/volumeunassigned',K.ADAPTER,R.STRUCTURAL),stock('kit_preparation','SourceactualopaqueOTHERpre-PCRmolecularpreparationexcludingseparatelytrackedadapteraliquot; exactenzymes/amountsunassembled',K.MOLECULAR,R.STRUCTURAL)))
    prep=g.input('kit_pcr_preparation','SourceONEhistoricalkitPCRpreparation;internalstocks/quantitiesunassembled',K.MOLECULAR);poly=g.product('kit_polymerase','PolymerasecomponentofsamepreparedPCRmixture',K.CATALYST);nts=g.product('kit_nucleotides','MolecularnucleotidecomponentofsamepreparedPCRmixture;exactpoolunassigned',K.MOLECULAR);f=g.product('kit_forward','ForwardprimercomponentofsamePCRpreparation;sourceorderingunassembled',K.OLIGO);r=g.product('kit_reverse','ReverseprimercomponentofsamePCRpreparation;sourceorderingunassembled',K.OLIGO)
    g.react('pcr_components','RepresentcopresentcomponentsofONEsourcePCRpreparation,nophysicalfractionation/extrastockaddition',S.GENERIC,((prep,R.TEMPLATE),),(poly,nts,f,r),changes_molecule=False)
    op('pcr','SourceSEVENPCRcycles;sourcehistoricalprogrammeunassembled',S.PCR,((poly,R.CATALYST),(nts,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Actualsource7cyclesknown;historicalprimers/wholeproductsunassembled,notconfirmedmissing.'))
    op('kit_complete','Completehistoricalkitpost-PCRprocessing;internalworkupandfractionsunassembled',S.MAGIC,(cond('kit_postprocessing','Sourcehistoricalkitcompletion'),))
    qc('final_qc','SourceQubit/TapeStationanalyticalaliquots;separatefromproduction')
    aliquot('normalized_aliquot','Sourcenormalizedaliquotforequimolarlibrarypool;volumeunassigned','Untransferredownlibrarypreparation;amount/laterfateunassigned')
    other=g.input('other_library_stocks','Otherindependentlysource-preparedlibrariesforequimolarsequencingpool;preparation/provenanceoutsidecurrenttechnical-replicategraph',K.LIBRARY);oa=g.product('other_normalized_aliquots','Sourceothernormalizedlibraryaliquots',K.LIBRARY);ore=g.product('other_library_remainder','Untransferredotherlibrarystocks',K.LIBRARY);g.react('other_library_aliquots','Sourceactualother-libraryaliquotselectionforpool',S.SPLIT,((other,R.TEMPLATE),),(oa,ore));g.mark_terminal(Terminal(ore.id,Disposition.UNRESOLVED,'Otherlibrarystockremainder,notsilentlyconsumed.'))
    lib=g.product('library','SourceequimolarpooledNovaSeqsequencinglibrary',K.LIBRARY);g.react('pool','Actualequimolarpoolofselectedownandothersourcealiquots',S.POOL,((x,R.TEMPLATE),(oa,R.TEMPLATE)),(lib,))
    g.sequence(lib,SequencingPlan(limitation='PublishedsourceNovaSeqPE50. Biologicalreplicateincludespooledtechnicaltriplicates;pooltiming/technicalreplicatesindividualcontributionsunassigned,nofabricatedearlypool. Wholeadapter/index/readlandingunassembled. Coverage/DPCratio/readfiltersarecomputational,notphysicalDNAremoval.'))
    g.add_source_item(SourceItem('pending','CitedProtocolExchangev1accessactuallyCloudflareblocked;historicalkitoligos/endsandtechnicalreplicatepooltimingunassembled. Table6CRISPRi/otherassayoligosnotimposedonDPCkit.','Selectedsource-positiveprivatepreparation;noabsence/fullsourceaudit.'));return g
def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Inputcontrolisrealextractedmaterial',g.terminals['input_dna'].disposition,Disposition.UNRESOLVED)
    c('Solublequantificationdoesnotconsumeentirestock',g.terminals['soluble_remainder'].disposition,Disposition.UNRESOLVED)
    c('Threeactualreprecipitationwashes',[g.steps['wash'+str(i)].inputs[0].material for i in(1,2,3)],['dpc_pellet','wash1_kept','wash2_kept'])
    c('BSAaddedafterProteinaseK',g.steps['bsa'].inputs[0].material,'proteinase_product')
    c('CollectedDPCDNAsupernatantbeforeRNase',g.steps['rnase'].inputs[0].material,'released_supernatant')
    c('50ngaliquotfromquantifiedDPCfraction',g.steps['library_dna'].inputs[0].material,'qubit_kept')
    c('Firstglycogencarriercontributesmolecularmaterial',[(p.material,p.role)for p in g.steps['precipitate_dna1'].inputs if p.material=='glycogen1'],[('glycogen1',R.STRUCTURAL)])
    c('SYBRGoldcontributesmolecularmaterial',[(p.material,p.role)for p in g.steps['gel_stain'].inputs if p.material=='sybr_gold'],[('sybr_gold',R.STRUCTURAL)])
    c('SpinXgelresidualexplicit',g.terminals['spinx_residual'].disposition,Disposition.UNRESOLVED)
    c('PCRconsumesproducedsame-stockcomponents',[p.material for p in g.steps['pcr'].inputs if p.role==R.STRUCTURAL],['kit_nucleotides','kit_forward','kit_reverse'])
    c('Poolusesnormalizedownaliquot',g.steps['pool'].inputs[0].material,'normalized_aliquot_used')
    c('Poolusesexplicitotheraliquots',g.steps['pool'].inputs[1].material,'other_normalized_aliquots')
    c('Untransferredownlibraryisretainedasremainder',g.terminals['normalized_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    c('Sourceinitialfragmentationonly',sum(s.kind==S.FRAGMENTATION for s in g.steps.values()),1)
    c('OtherprePCRpreparationisactualmolecularcontribution',[(p.material,p.role)for p in g.steps['kit_pre_pcr'].inputs if p.material=='kit_preparation'],[('kit_preparation',R.STRUCTURAL)])
    c('Privatepartialscope',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
