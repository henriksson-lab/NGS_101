"""Source-local TF-chRDP: shared immunocapture, independent DNA/RNA libraries."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceCoverage,SequencingPlan
from pcr import PCRGeometry
AUTHORITY='Published2025PMC11925928 main/reagents/Methods265–479/Results497–503 andactualFigure1PDFp9. ONEoriginal; no vendor/predecessor/wholekitoligo authority imported. HCT116/5FU/p53 selectedreference route.'

PAGE_FILE = 'tf-chrdp.html'

# Curator-reviewed inventory of this source's actual inputs and operations.
# These fixed source groups do not enumerate graph nodes at runtime. Repetition
# follows the source's control/production arms, not a shared biological recipe.
SHARED_CAPTURE_SCOPE = (
    'source_cells','treat','fu','initial_wash','initial_wash_conditions',
    'fa','fa_stock','fa_quench','fa_glycine','fa_quench_conditions',
    'fa_wash','fa_wash_conditions','egs','egs_stock','egs_quench',
    'egs_glycine','egs_quench_conditions','egs_wash','egs_wash_conditions',
    'treatment_pool','treatment_pool_scope','sample','lysis','lysis_buffer',
    'permeabilize','sds','nuclear_wash','nuclear_wash_conditions',
    'sonic_resuspend','sonic_buffer','sonic_allocate','sonic1','sonic1_conditions',
    'sonic2','sonic2_conditions','sonic3','sonic3_conditions','sonic_pool',
    'clarify','clarify_conditions','preclear_beads_raw','preclear_beads_wash',
    'preclear_beads_wash_conditions','ip_beads_raw','ip_beads_wash',
    'ip_beads_wash_conditions','anti_p53','antibody_load','antibody_conditions',
    'preclear','preclear_remove','preclear_remove_conditions','preip_qc_allocate',
    'antibody_wash','antibody_wash_conditions','ip','ip_wash','ip_wash_conditions',
    'three_fractions',
)
DNA_RECOVERY_SCOPE = ('decrosslink','pk','buffers','column','column_conditions','collect','quant')
RNA_RECOVERY_SCOPE = (
    'decrosslink','pk','rna_buffer','boil','boil_conditions','trizol','trizol_supply',
    'chloroform','chloroform_supply','chloroform_conditions','ethanol','ethanol_supply',
    'column_load','column_load_conditions','dnase_mix_stock','dnase_mix_components',
    'dnase','rwt_first','rwt_conditions','rwt_reapply','rwt_reapply_conditions',
    'rna_wash','rna_wash_conditions','dnase_inactivate','dnase_stop','rt_input',
    'hexamers','rt_dntp','prime','ssiii','rnaseout','rt_other','rt','pol1','rnaseh',
    'ligase','second_dntp','second_buffer','second','t4pol','t4pol_stock',
    'edta','stop','stop_conditions','cdna_clean','cdna_clean_conditions',
    'cdna_collect','cdna_qc',
)
TAGMENT_PCR_SCOPE = (
    'template','tn5_stock','tn5_components','tag','tag_conditions','sds','sds_supply',
    'master_stock','master_components','index_one','index_two','pcr',
)
TESTING_SCOPE = ('test_clean','test_clean_conditions','test_analysis')
LIBRARY_SELECTION_SCOPE = (
    'volume','eb','first_size','first_selection','second_size','second_size_conditions',
    'second_selection','library_collect','library_quant','loading','load','library',
)

def declare_page_scope(g):
    g.declare_reviewed_scope(SHARED_CAPTURE_SCOPE,
        'Published TF-chRDP main Methods3.1–3.9: selected HCT116/5-FU/p53 route, same-treatment pool, actual three sonication tubes, pre-IP QC aliquots, bead preparation and three post-IP fractions. Aggregate cardinality and molecular occupancy remain unassigned.')
    for prefix in ('dna_control','dna'):
        g.declare_reviewed_scope(tuple(prefix+'_'+target for target in DNA_RECOVERY_SCOPE),
            'Published Methods3.11.1/3.12: separately processed DNA QC and IP-DNA, source buffers/proteinase K, column workup, eluate transfer and analytical fractions.')
    g.declare_reviewed_scope(('dna_bead_remove','dna_bead_remove_conditions'),
        'Published Methods3.12: magnetic removal of production IP beads before DNA purification; no bead-removal step invented for the bead-free input control.')
    for prefix in ('rna_control','rna'):
        g.declare_reviewed_scope(tuple(prefix+'_'+target for target in RNA_RECOVERY_SCOPE),
            'Published Methods3.11.2/3.14–3.17: distinct RNA workup, actual retained/reapplied RWT flow-through, DNase mix, random-hexamer SSIII and second-strand/T4-polymerase processing. Inventory classes do not assign full kit formulations or occupied products.')
    for prefix in ('dna_control','dna','rna_control','rna'):
        g.declare_reviewed_scope(tuple(prefix+'_'+target for target in TAGMENT_PCR_SCOPE),
            'Published Methods3.13.1–3.13.2 and explicit RNA reuse3.18: actual Tn5/PCR supplies, same-stock component classes and paired index-primer functions. DNA-test dose is a reference, not an invented optimized production amount or canonical oligo body.')
    for prefix in ('dna_control','rna_control'):
        g.declare_reviewed_scope(tuple(prefix+'_'+target for target in TESTING_SCOPE),
            'Published Methods3.13.1/3.18: input-derived test cleanup and Qubit/Qsep suitability measurement; controls are not sequenced or physically pooled with production libraries.')
    for prefix in ('dna','rna'):
        g.declare_reviewed_scope(tuple(prefix+'_'+target for target in LIBRARY_SELECTION_SCOPE),
            'Published Methods3.13.3–3.13.4 and3.18: first0.61x supernatant, second19ul bead addition, eluate/analytical/loading allocations and independent NovaSeqPE150 DNA/RNA endpoints. Optimized yields, whole molecules and read landing remain unresolved.')
    for item,target in (('scope','source_cells'),('array_scope','sonic_allocate'),('kit_scope','dna_tag')):
        g.cover_source_item(SourceCoverage(item,'modelled',target))
    return g.require_publishable()

def unknown_scene(name):
    sc=Scene();sc.strand('source',[Segment(name,'[whole source molecule unassembled]',placeholder=True)]);return sc

class T:
    def __init__(self,g,x):self.g=g;self.x=x
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return self.g.input(key,label,kind,**kw),role
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,ports=(),kind=S.MAGIC,geometry=None):
        y=self.g.product(key+'_product','Actual '+label+' mixture; wholeoccupiedmolecules/recoveryunassembled',K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Source-localpositiveoperation; completeoccupiedsubstrates/products/formulations/individualmolecularfatesunassembled, notabsence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
    def clean(self,key,label,ports=()):
        y=self.g.product(key+'_kept','Actual retained '+label+' fraction; recoveryunassigned',K.MOLECULAR);w=self.g.product(key+'_excluded','Actual removed '+label+' fractions',K.WASTE);self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_conditions',label))+tuple(ports),(y,w));self.g.discard(w,'Actualsourceexcludedfraction, notperfectmolecularspecificity');self.x=y;return y
    def take(self,key,label,rest,qc=False):
        a=self.g.product(key+'_selected',label,K.QC if qc else K.MOLECULAR);r=self.g.product(key+'_rest',rest,K.MOLECULAR);self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r));self.g.mark_terminal(Terminal(a.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else rest));self.x=r if qc else a;return self.x
    def components(self,key,label,items):
        st=self.g.input(key+'_stock',label,K.MOLECULAR);xs=[self.g.product(key+'_'+name,caption,kind,**kw)for name,caption,kind,role,kw in items];self.g.react(key+'_components','Components ofSAMEactualstock; nofreshduplicateadditions/physicalseparation',S.GENERIC,((st,R.TEMPLATE),),tuple(xs),changes_molecule=False);return tuple((x,row[3])for x,row in zip(xs,items))

def washed_beads(g,key,label):
    raw=g.input(key+'_raw',label,K.MOLECULAR);b=T(g,raw);b.clean(key+'_wash','Three1ml0.1percentTritonDPBS washes; retainProteinGbeads');return b.x

def dna_recover(t,key,control=False):
    t.op(key+'_decrosslink','65Covernight900rpm DNAdecrosslink',ports=(t.stock(key+'_pk','Actual10ulProteinaseK'),t.stock(key+'_buffers','Actual'+('30ulTE+50ulChIPelutionbuffer into10ulQC'if control else'40ulTE+50ulChIPelutionbuffer intobeads')+'; actualsourceallocations',K.MOLECULAR,R.STRUCTURAL)))
    if not control:t.clean(key+'_bead_remove','Magnetic2min retainDNAsupernatant, excludeoriginalIPbeads')
    t.clean(key+'_column','ZymoGenomicDNACleanConcentrator10:fivevolbindingbuffer/reapplyflowthrough/two200ulwashes/elute15ulEB')
    t.take(key+'_collect','Actual14ulof15ulDNAeluate','Untransferred1ulDNAelutionstock')
    t.take(key+'_quant','ActualQubit4 analyticalallocation; QsepfragmentQC oncontrol only','DNAstockafteractualanalyticalallocation',qc=True)
    return t.x

def rna_to_cdna(t,key,control=False):
    t.op(key+'_decrosslink','50C45min900rpm RNAdecrosslink',ports=(t.stock(key+'_pk','Actual5ulProteinaseK'),t.stock(key+'_rna_buffer','Actual'+('85ulRNAProteinaseKbuffer into10ulRNAQC'if control else'95ulRNAProteinaseKbuffer intoRNA-IPbeads')+' sourcebufferTris/NaCl/EDTA/SDS',K.MOLECULAR,R.STRUCTURAL)))
    t.op(key+'_boil','95C10min thenice5min',ports=(t.cond(key+'_boil_conditions','ActualsourceRNAroute, notDNA65Covernight'),))
    t.op(key+'_trizol','Add500ulTRIzol/vortex10s/RT10min',ports=(t.stock(key+'_trizol_supply','Actual500ulTRIzol',K.MOLECULAR,R.STRUCTURAL),))
    t.clean(key+'_chloroform','Add100ulchloroform/16100rcf4C15min; RETAIN400ulaqueous fraction',ports=(t.stock(key+'_chloroform_supply','Actual100ulchloroform allocation',K.MOLECULAR,R.STRUCTURAL),))
    t.op(key+'_ethanol','Actual600ul100percentEtOH',ports=(t.stock(key+'_ethanol_supply','Actual600ulEtOH',K.MOLECULAR,R.STRUCTURAL),))
    t.clean(key+'_column_load','miRNeasymicro sequential≤700ulcolumnloading and350ulRWTwash; retainRNA-boundcolumn')
    dn=t.components(key+'_dnase_mix','Actual80ulpreparedDNaseIincubationmix:70ulRDD+10ulDNasestock',(('dnase','DNaseIcomponentofsame80ulphysicalmix',K.CATALYST,R.CATALYST,{}),('buf','RDDcomponentofsame80ulphysicalmix',K.CONDITION,R.CONDITION,{})))
    t.op(key+'_dnase','On-columnDNase20–30C15min',ports=dn,kind=S.DIGESTION)
    # Source specifically reapplies the RWT flow-through once before discarding it.
    bound=t.g.product(key+'_rwt_bound','Actualcolumn-boundRNAfraction afterfirst500ulRWTpass',K.MOLECULAR);flow=t.g.product(key+'_rwt_flow','ActualRWTflowthrough RETAINEDforone reapplication',K.MOLECULAR)
    t.g.react(key+'_rwt_first','Sourcefirst500ulRWTpass:keepflowthrough temporarily',S.CLEANUP,((t.x,R.TEMPLATE),t.cond(key+'_rwt_conditions','Source500ulRWT/10000rpm15s')),(bound,flow))
    rejoined=t.g.product(key+'_rwt_reapplied','ActualcolumnRNA aftersameflowthroughreapplied',K.MOLECULAR);w=t.g.product(key+'_rwt_discard','Actualexcludedsecondpassflowthrough',K.WASTE);t.g.react(key+'_rwt_reapply','Reapply SAMEretainedRWTflowthrough to SAMEcolumn',S.CLEANUP,((bound,R.TEMPLATE),(flow,R.STRUCTURAL),t.cond(key+'_rwt_reapply_conditions','10000rpm15s; secondflowthroughdiscarded')),(rejoined,w));t.g.discard(w,'Source secondpassRWTflowthroughdiscarded');t.x=rejoined
    t.clean(key+'_rna_wash','500ulRPE then500ul80percentEtOH/2minspin/dry5min/elute14ulwater')
    t.op(key+'_dnase_inactivate','65C15min afterRNAelution',ports=(t.cond(key+'_dnase_stop','Actualsourceheatprogramme'),),kind=S.QUENCHING)
    t.take(key+'_rt_input','Actual13ulRNAof14uleluate','Untransferred1ulRNAeluate')
    n6=Scene();n6.strand('randomprimer',[Segment('Randomhexamer','NNNNNN',note='Source randomprimingfamily; NOT deduplicationUMI')])
    t.op(key+'_prime','Add1ul50ngulrandomhexamers+1ul10mMdNTPs;65C5min/ice≥1min',ports=(t.stock(key+'_hexamers','Actual1ulrandomhexamerstock; sugar/endmods/wholefamilyunassembled',K.OLIGO,R.STRUCTURAL,state=n6),t.stock(key+'_rt_dntp','Actual1ul10mMdNTPmix',K.MOLECULAR,R.STRUCTURAL)))
    t.op(key+'_rt','Actual30ulfirststrandreaction25C10min/50C50min',ports=(t.stock(key+'_ssiii','Actual1ul200UulSuperScriptIII'),t.stock(key+'_rnaseout','Actual1ul40UulRNaseOUTprotein inhibitor, not addednucleotides',K.MOLECULAR,R.STRUCTURAL),t.stock(key+'_rt_other','Actual4ulwater/3ul10xRTbuffer/4ul25mMMgCl2/2ul0.1MDTT supplies; notclaimedpremix',K.MOLECULAR,R.STRUCTURAL)))
    t.op(key+'_second','Source150ulsecondstrand16C2h',ports=(t.stock(key+'_pol1','Actual4ul10UulEcoliDNApolymeraseI'),t.stock(key+'_rnaseh','Actual1ul2UulEcoliRNaseH'),t.stock(key+'_ligase','Actual1ul10UulEcoliDNAligase'),t.stock(key+'_second_dntp','Actual3ul10mMdNTPmix',K.MOLECULAR,R.STRUCTURAL),t.cond(key+'_second_buffer','Actual30ulfirststrand+81ulwater+30ul5xSecondStrandbuffer; temperature≤16C')))
    t.op(key+'_t4pol','Add3.3ul3UulT4DNApolymerase12C15min1100rpm',ports=(t.stock(key+'_t4pol_stock','Actual3.3ulT4DNApolymerase'),))
    t.op(key+'_stop','Add10ul0.5MEDTAonice',ports=(t.stock(key+'_edta','Actual10ul0.5MEDTA',K.MOLECULAR,R.STRUCTURAL),t.cond(key+'_stop_conditions','Actualsourceiceprogramme')),kind=S.QUENCHING)
    t.clean(key+'_cdna_clean','Actual294ulAMPure(1.8x)/two800ul80percentEtOH/dry~3min/elute14ulEB')
    t.take(key+'_cdna_collect','Actual12ulof14ulcDNAeluate','Untransferred2ulcDNAelutionstock')
    t.take(key+'_cdna_qc','ActualQubit4/Qsep analyticalallocation','cDNAstockafteractualQC',qc=True)
    return t.x

def tagment(t,key,testing=False):
    dose='ActualspecifiedDNAQCtest8ul4ngtemplate/0.33ulTn5'if key=='dna_control'else'Source-selectedworkingtemplate/Tn5allocation; actualoptimizedamountunassigned. SourceDNA-testingreference8ul4ng/.33ulTn5 isnotassignedasoptimizeddose'
    t.take(key+'_template',dose,'UntransferredsourceDNA/cDNAstock afterselectedTn5allocation')
    loaded=t.components(key+'_tn5','ActualsourceTn5stockBGILS-EZ-E-00009O selectedallocation; '+dose+'; wholeloadedgeometryunassembled',(('enzyme','Tn5catalystcomponentofSAMEactualselectedtransposomestock; sourceDNA-testreference0.33ul, optimizeddoseunassigned',K.CATALYST,R.CATALYST,{}),('adapters','SourcefunctionalDNAadaptercomponentsofSAMEactualselectedTn5stock; wholeorders/speciesunassembled',K.ADAPTER,R.STRUCTURAL,{'state':unknown_scene('SourceTn5adapterinventoryunassembled')})))
    t.op(key+'_tag','Source55C10min/4C Tn5tagmentation',ports=loaded+(t.cond(key+'_tag_conditions','Actual6.67ulwater+8ultemplate+5ul4xTHSTD+0.33ulTn5reference20ul; cDNArecipe reuseexplicit3.18, actualoptimizedinputconventionunassigned'),))
    t.op(key+'_sds','Add3ul0.2percentSDS/37C5min/RT5min/ice>5min',ports=(t.stock(key+'_sds_supply','Actual3ul0.2percentSDS',K.MOLECULAR,R.STRUCTURAL),))
    master=t.components(key+'_master','Actual25ulHigh-Fidelity2xPCRMasterMixNEBM0541 referenceallocation',(('enzyme','DNApolymerasecomponentofSAMEactual25ulmastermix',K.CATALYST,R.CATALYST,{}),('nt','NucleotidecomponentsofSAMEactual25ulmastermix',K.MOLECULAR,R.STRUCTURAL,{}),('buffer','BuffercomponentsofSAMEactual25ulmastermix',K.CONDITION,R.CONDITION,{})))
    ps=[]
    for number,role in [('one','i7'),('two','i5')]:
        sc=Scene();sc.strand('primer',[Segment('SourceIndex'+number,'[whole indexed primer unassembled]',placeholder=True,feature=feature('tfchrdp.sample.'+role,'sample_index','unknown',note='Source TruePrepV2Index'+number+' '+role+'; selectedvalue/width/wholeorderunassembled'))]);ps.append(t.stock(key+'_index_'+number,'Actual1ulIndex'+number+' '+role+' referenceprimerstock; wholeordersunassembled',K.OLIGO,R.STRUCTURAL,state=sc))
    t.op(key+'_pcr','72C3min/98C30s/eightx98C15s60C30s72C3min/72C5min',ports=tuple(ps)+master,kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene('Whole taggedDNA/cDNAunassembled'),reason='KnownsourcePCRpair/classes/programme; wholevendorindex/adapters/occupiedsites/products unassembled. No canonicalsequence orcommonDNA/cDNAproductinherited.'))
    if testing:
        t.clean(key+'_test_clean','Source1xAMPureaftertestingPCR')
        t.take(key+'_test_analysis','ActualQubit/QseptestingQC; source250–600bp/>30ngsuitabilitycriteria, not imposed yield','UntransferredTn5teststockafteranalysis',qc=True);t.g.mark_terminal(Terminal(t.x.id,Disposition.UNRESOLVED,'Remaininganalyticalteststock, notsequenced'));return
    t.op(key+'_volume','EBadjustPCRproductto100ul',ports=(t.stock(key+'_eb','ActualEBadjustmentallocation; amountunassigned',K.MOLECULAR,R.STRUCTURAL),),kind=S.GENERIC)
    sup=t.g.product(key+'_first_supernatant','ActualRETAINEDfirst0.61xsupernatant',K.MOLECULAR);beads=t.g.product(key+'_first_beads','Actualexcludedfirstbeadboundfraction',K.WASTE);t.g.react(key+'_first_size','61ulAMPureto100ulPCR/10minrotation/5minmagnet; RETAINsupernatant',S.CLEANUP,((t.x,R.TEMPLATE),t.stock(key+'_first_selection','Actual61ulAMPure',K.SELECTION,R.SELECTION)),(sup,beads));t.g.discard(beads,'Actualfirstbeadboundfractionexcluded');t.x=sup
    t.clean(key+'_second_size','Actual19ulAMPure(.8-.61)xORIGINAL100ul; RETAINbeads/two800ul80percentEtOH/elute8ulEB',ports=(t.stock(key+'_second_selection','Actual19ulAMPure',K.SELECTION,R.SELECTION),))
    t.take(key+'_library_collect','Actual6ulof8ulfinishedlibraryeluate','Untransferred2ulfinishedlibraryeluate')
    t.take(key+'_library_quant','ActualQubit4analyticalallocation','Finallibraryremainingafterquantification',qc=True)
    t.take(key+'_loading','Selectedfinalsequencingloadingaliquot; amountunassigned','UntransferredQCedfinishedlibrarystock')
    lib=t.g.product(key+'_library','Actualsource'+key+' libraryloadingmixture; wholeconstructs/readlandingunassembled',K.LIBRARY);t.g.react(key+'_load','SourceNovaSeq6000PE150',S.GENERIC,((t.x,R.TEMPLATE),),(lib,),changes_molecule=False);t.g.sequence(lib,SequencingPlan(read_lengths={'R1':150,'R2':150},limitation='SourceDNA/RNAseparatePE150 libraries; wholeTn5/indexorders/occupiedgenomicorRNA-copyproducts/landingunassembled. No physicalDNA/RNApoolinvented.'))

def protocol_graph():
    g=ProtocolGraph('TF-chRDP sharedp53capture and separateDNA/RNAlibraryroutes',authority=AUTHORITY,representation='published_instance');cells=g.input('source_cells','Actualsame-treatmentHCT116culturecollection; dishnumbersunassigned, notone30Mcelldish',K.BIOLOGICAL,primary=True);t=T(g,cells)
    t.op('treat','375uM5FU9h selectedtreatment, notDMSOcontrol',ports=(t.stock('fu','Actual5FU375uM culturemediumallocation; aggregateamountunassigned',K.MOLECULAR,R.STRUCTURAL),))
    t.clean('initial_wash','Two20mlDPBS washes perdishreference')
    t.op('fa','FAFIRST:1percent20mlperdishreference/20minRT50rpm',ports=(t.stock('fa_stock','ActualaggregateFAstockallocations;20ml1percentperdishreference, multiplicityunassigned',K.MOLECULAR,R.STRUCTURAL),))
    t.op('fa_quench','Glycine1.74ml2.5Mperdishreference/0.2Mfinal/10min',ports=(t.stock('fa_glycine','Actualaggregateglycineallocations; sourceperdishreference retained',K.MOLECULAR,R.STRUCTURAL),t.cond('fa_quench_conditions','ActualRT10min50rpm')),kind=S.QUENCHING);t.clean('fa_wash','Two20mlDPBS washesperdishreference')
    t.op('egs','EGSSECOND:2mM20mlperdishreference/45minRT50rpm',ports=(t.stock('egs_stock','ActualaggregateEGSstockallocations, dishnumbersunassigned',K.MOLECULAR,R.STRUCTURAL),))
    t.op('egs_quench','Glycine1.74ml2.5Mperdishreference/10min',ports=(t.stock('egs_glycine','Actualaggregateglycineallocations',K.MOLECULAR,R.STRUCTURAL),t.cond('egs_quench_conditions','ActualRT10min50rpm')),kind=S.QUENCHING);t.clean('egs_wash','TwoDPBS washes then coldDPBS/PI rinse/scrape')
    t.op('treatment_pool','Physicalsame-treatmentcellpool after~20Mcellsperdishcollection; dishcardinalityunassembled',ports=(t.cond('treatment_pool_scope','Actualsourcephysicalpool boundary; unknownnumberofcontributingdishstocks, nofakeone-dish30Msample'),))
    t.take('sample','Actual~30Mcellselectedsample;5000rpm5minpellet/freezeoptional','Untransferredsame-treatmentcellpool')
    t.op('lysis','10ml0.1percentSDSFAbuffer+PI/RT15min',ports=(t.stock('lysis_buffer','Actual10mlSDSFAbufferwithPI allocation',K.MOLECULAR,R.STRUCTURAL),))
    t.op('permeabilize','Add1ml10percentSDS/37C20min',ports=(t.stock('sds','Actual1ml10percentSDS',K.MOLECULAR,R.STRUCTURAL),))
    t.clean('nuclear_wash','2000rpm10min4C,retainpellet;10ml0.1percentSDSFAbufferwash then12000rpm15min')
    t.op('sonic_resuspend','Resuspendpellet3mlSDSlysisbuffer/PI',ports=(t.stock('sonic_buffer','Actual3ml0.1percentSDSlysisbuffer/PI',K.MOLECULAR,R.STRUCTURAL),))
    tubes=[g.product('sonic_tube_'+str(i),'Actual1mlsonicationtube'+str(i),K.MOLECULAR)for i in range(1,4)];g.react('sonic_allocate','Actual3mlstockdividedintoTHREE1mlsonicationtubes',S.SPLIT,((t.x,R.TEMPLATE),),tuple(tubes));sonics=[]
    for i,x in enumerate(tubes,1):
        t.x=x;t.op('sonic'+str(i),'Source38percentamplitude5min/30son-off/onice',ports=(t.cond('sonic'+str(i)+'_conditions','Actualsameprogramme; fragment/enddistributionunassigned'),),kind=S.FRAGMENTATION);sonics.append(t.x)
    pool=g.product('sonicated_pool','Actualsonicatesreturnedtooriginal15mltube',K.MOLECULAR);g.react('sonic_pool','Sourcephysicalrecombination AFTERthreeactualsonications',S.POOL,tuple((x,R.TEMPLATE)for x in sonics),(pool,));t.x=pool;t.clean('clarify','2000rpm15min4C,RETAINchromatinsupernatant; source1mlcollectiontubes, exactrecoveredarraydistributionunassembled')
    pre=washed_beads(g,'preclear_beads','ActualaggregateProteinGpreclearbeadstock;50ulper-samplereference, subsequent1mlchromatinarray conventionunassembled');cap=washed_beads(g,'ip_beads','ActualaggregateProteinGIPbeadstock;150ulperantibodyreference')
    ab=g.input('anti_p53','Actual10uganti-p53sc126persamplereference',K.MOLECULAR);loaded=g.product('antibody_beads','Actualanti-p53ProteinGbeadpreparation after8h4C12rpm/500ulTritonDPBSreference',K.MOLECULAR);g.react('antibody_load','Bindanti-p53toSAMEprewashedIPbeads',S.MAGIC,((cap,R.TEMPLATE),(ab,R.STRUCTURAL),t.cond('antibody_conditions','Source8h4C; actualaggregatebead/antibodydose-arrayconventionunassembled')),(loaded,),reason='Actualsourceantibody-loadedselectionstock; fullprotein/crosslinkoccupancy unassembled.')
    t.op('preclear','Actualprewashedbeads+chromatin4C≥1h12rpm',ports=((pre,R.STRUCTURAL),))
    t.clean('preclear_remove','Retainpreclearedchromatin, remove nonspecificProteinGboundfraction')
    dq=g.product('dna_qc','Actual10ulpreIPDNAQCsample',K.MOLECULAR);rq=g.product('rna_qc','Actual10ulpreIPRNAQCsample',K.MOLECULAR);production=g.product('ip_chromatin','Actualpreclearedchromatinremainingaftertwo10ulQCaliquots',K.MOLECULAR);g.react('preip_qc_allocate','SourceTWO10ulcontrols BEFOREp53IP',S.SPLIT,((t.x,R.TEMPLATE),),(dq,rq,production))
    before=t.x;t.x=loaded;t.clean('antibody_wash','Two1ml0.1percentTritonDPBS washesbeforechromatinIP');washed=t.x;t.x=production
    t.op('ip','Actualp53IPovernight4C12rpm',ports=((washed,R.STRUCTURAL),))
    t.clean('ip_wash','RemoveIPsupernatant, threeLowSalt/twoHighSalt/twoTEwashes1mleach4C5min')
    dna=g.product('dna_fraction','ActualpostwashbeadcomplexallocationforDNA',K.MOLECULAR);rna=g.product('rna_fraction','ActualpostwashbeadcomplexallocationforRNA',K.MOLECULAR);protein=g.product('protein_fraction','ActualpostwashbeadcomplexallocationforLC-MS/MSproteinidentification; NOTDNAbarcode',K.MOLECULAR);g.react('three_fractions','AllocateactualwashedcomplexestoTHREEseparatetubes; fractionsizesunassigned',S.SPLIT,((t.x,R.TEMPLATE),),(dna,rna,protein));g.mark_terminal(Terminal(protein.id,Disposition.UNRESOLVED,'Actualnon-NGSassaycontinuation outside thismodel; LC-MS/MSfraction isnot aQCcontrol orDNA/RNAsequencinglibrary'))
    for key,stock in [('dna_control',dq),('dna',dna)]:
        t.x=stock;dna_recover(t,key,control=key.endswith('control'));tagment(t,key,testing=key.endswith('control'))
    for key,stock in [('rna_control',rq),('rna',rna)]:
        t.x=stock;rna_to_cdna(t,key,control=key.endswith('control'));tagment(t,key,testing=key.endswith('control'))
    g.add_source_item(SourceItem('scope',AUTHORITY,'Source-qualified shared IP and independent DNA/RNA workflows; no perfect-purity or whole-library molecule claim.'))
    g.add_source_item(SourceItem('array_scope','~20Mcellsperdish pooledbytreatment→30Mselectedsample; threeactual1mlsonications.','Actualdishcardinality, recovered1mlchromatindistribution andper-samplebeadstockmultiplicityunassembled; nofakephysicalcell/nucleic-acidpool orone-dosecloning.'))
    g.add_source_item(SourceItem('kit_scope','SourceTn5DNAadapterinstallation and25ulMasterMixM0541 arepositive; actualwholekitoligosunassembled.','Knownfunctionaladapter/enzyme andPCRenzyme/nucleotide/bufferclassesfromSAMEactualstocks; not freshindependentadditions or borrowedcanonicalends.'));return declare_page_scope(g)

def render_page():
    from page import head,info,caveat,table
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph()
    branches=table(('Fraction','Preparation and endpoint'),(
        ('Captured DNA','Proteinase K / overnight at 65°C workup, DNA purification, source-calibrated Tn5/PCR and two-stage bead selection; NovaSeq DNA library.'),
        ('Captured RNA','Milder 50°C proteinase treatment, heat/TRIzol and column workup, DNase, random-hexamer RT and second-strand cDNA synthesis; separate Tn5/PCR RNA-derived library.'),
        ('Captured protein','A separate LC-MS/MS allocation. Its detailed proteomics continuation is outside this DNA/RNA model.'),
    ))
    reads=table(('Library','R1 cycles','R2 cycles'),[(name.replace('_library','').upper(),str(plan.read_lengths['R1']),str(plan.read_lengths['R2']))for name,plan in g.sequencing.items()])
    return assemble_page(g,[head('TF-chRDP — joint capture, separate DNA/RNA libraries'),
        '<main class="wrap" data-protocol-complete="true"><h1>TF-chRDP</h1>',
        info('<a href="https://doi.org/10.3389/fcell.2025.1561540">Ning et al., Frontiers in Cell and Developmental Biology (2025)</a> · <a href="01_tf-chrdp.html">Research notes</a>'),
        '<p>Capture transcription-factor-associated chromatin complexes once, then divide them into DNA, RNA and protein fractions. This page follows the published HCT116/5-FU/p53 preparation.</p>',
        caveat('The source-supported preparation and physical fractions are represented. Complete substrates, crosslinks, loaded adapters, index-primer sequences and occupied read sites remain unassembled. The optimized production Tn5 amount is unassigned; the DNA-control test dose is kept as a reference.'),
        '<h2>Shared capture</h2><ol><li>Treat cells with 5-FU, then crosslink with formaldehyde followed by EGS, with separate glycine quenches and washes.</li><li>Pool same-treatment cells, select a 30-million-cell sample, lyse/permeabilize nuclei and divide the sonication suspension into three 1 ml tubes. Recombine after sonication.</li><li>Preclear chromatin and retain separate 10 µl DNA and RNA controls. Immunoprecipitate the remaining chromatin with prepared anti-p53/Protein G beads.</li><li>Wash the captured complexes and allocate separate DNA, RNA and protein fractions.</li></ol>',
        '<h2>Separate downstream preparations</h2>',branches,
        '<h2>Library construction</h2><p>The source explicitly reuses its DNA tagmentation/PCR procedure for RNA-derived double-stranded cDNA. The same actual transposome inventory supplies catalyst and adapter classes; the actual PCR master mix supplies enzyme, nucleotide and buffer classes. Whole adapter/index orders are not filled from a current kit or platform name.</p>',
        '<p>Production libraries retain the first 0.61× bead-selection supernatant, then add 19 µl more beads relative to the original 100 µl preparation and retain the bead-bound fraction. Input-derived test reactions and analytical aliquots stay separate from the sequencing endpoints.</p>',
        '<h2>Read setup</h2><p>Illumina NovaSeq 6000, paired-end 150 bases, for each of the independent DNA and RNA-derived libraries. The model does not assert a physical pool between them.</p>',reads,
        PlanLimitationSection('dna_library'),PlanLimitationSection('rna_library'),
        '<details><summary>Detailed workflow, controls and unresolved chemistry</summary>',WorkflowSection('Source-qualified capture and library workflows'),'</details></main>'])


def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('FA precedes EGS via actual washed stock',g.steps['egs'].inputs[0].material,'fa_wash_kept')
    c('Three actual sonication aliquots',g.steps['sonic_allocate'].outputs,('sonic_tube_1','sonic_tube_2','sonic_tube_3'))
    c('Sonication products recombined onlyafter separate actualtreatment',[x.material for x in g.steps['sonic_pool'].inputs],['sonic'+str(i)+'_product'for i in range(1,4)])
    c('SharedpreIPQCtwo10ulactualcontrols',[x for x in g.steps['preip_qc_allocate'].outputs],['dna_qc','rna_qc','ip_chromatin'])
    c('IP receives QC-depletedchromatin not QCaliquot',g.steps['ip'].inputs[0].material,'ip_chromatin')
    c('ActualthreeIPfractionsnotindividualclonebarcodes',g.steps['three_fractions'].outputs,('dna_fraction','rna_fraction','protein_fraction'))
    c('Protein is nonNGSassaycontinuation not a QCcontrol',g.terminals['protein_fraction'].disposition,Disposition.UNRESOLVED)
    c('DNA begins actual DNAfraction',g.steps['dna_decrosslink'].inputs[0].material,'dna_fraction')
    c('RNA begins actual RNAfraction',g.steps['rna_decrosslink'].inputs[0].material,'rna_fraction')
    c('SeparateDNAQCprocessed from actual10ulcontrol',g.steps['dna_control_decrosslink'].inputs[0].material,'dna_qc')
    c('SeparateRNAQCprocessed from actual10ulcontrol',g.steps['rna_control_decrosslink'].inputs[0].material,'rna_qc')
    for k in ('rna','rna_control'):
        c('DNase fromsame actual80ulmix '+k,g.producer_of(k+'_dnase_mix_dnase').id,k+'_dnase_mix_components')
        c('Actualchloroform supply port '+k,any(x.material==k+'_chloroform_supply'for x in g.steps[k+'_chloroform'].inputs),True)
        c('RWTflowthrough physically reapplied not destroyedfirst '+k,g.steps[k+'_rwt_reapply'].inputs[1].material,k+'_rwt_flow')
        c('Firststrand uses selected13ulof14 '+k,g.steps[k+'_prime'].inputs[0].material,k+'_rt_input_selected')
        c('N6randomprimerhasnoUMIfeature '+k,g.materials[k+'_hexamers'].state.strands['randomprimer'].segs[0].feature,None)
        c('Secondstrand not fresh unrelatedRNA '+k,g.steps[k+'_second'].inputs[0].material,k+'_rt_product')
        c('T4pol follows actual secondstrand '+k,g.steps[k+'_t4pol'].inputs[0].material,k+'_second_product')
    for k in ('dna_control','dna','rna_control','rna'):
        c('ActualTn5enzyme andadapter share physicalstock '+k,g.producer_of(k+'_tn5_enzyme').id,g.producer_of(k+'_tn5_adapters').id)
        c('ActualPCRenzyme andNT share physicalmastermix '+k,g.producer_of(k+'_master_enzyme').id,g.producer_of(k+'_master_nt').id)
    for k in ('dna','rna'):
        c('Secondsizefraction uses retained .61xsupernatant '+k,g.steps[k+'_second_size'].inputs[0].material,k+'_first_supernatant')
        c('Libraryonlyselectedfinalloadingstock '+k,g.steps[k+'_load'].inputs[0].material,k+'_loading_selected')
        c('DNA test reference not fabricatedoptimizedproductiondose '+k,'actualoptimizedamountunassigned'in g.materials[k+'_template_selected'].label,True)
        c('Source PE150 '+k,g.sequencing[k+'_library'].read_lengths,{'R1':150,'R2':150})
    c('Onlytwoindependentproductionlibrariessequenced',set(g.sequencing),{'dna_library','rna_library'})
    c('DNA/RNA libraries not physicallymixed',any(set(x.material for x in st.inputs)=={'dna_library','rna_library'}for st in g.steps.values()),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'fcell.2025.1561540_PMC11925928.xml.txt'
    if have(c,Source(p,'https://doi.org/10.3389/fcell.2025.1561540')):
        t=p.read_text()
        for label,needle in [('ActualsourceQCpair','Take two 10 µL aliquots'),('Actualsource3fractions','aliquot the beads into 3 separate tubes'),('RNA source milderdecrosslink','50°C for 45 min'),('ActualsourceRWTreapplication','Reapply the flow-through to the RNeasy MinElute spin column'),('Actualsource19beads','19 μL ((0.8–0.61) × volume)'),('PositivePCRsamephysicalMasterMix','25 μL of High-Fidelity 2× PCR Master Mix')]:c(label,needle in t,True)
    c.report()
if __name__=='__main__':selftest()
