"""Private ChOR selected publishedEdU/qChOR capturetiming andphysicalstocks."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
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

AUTHORITY='PublishedPMC6202308main381–399 andpublishedSI1actualFigS1p2/FigS2p3. Two originals; notAMorlater2021recipe.'
def source_a_end():
    s=Scene();s.strand('top',[Segment('unknownbody','XXXX',placeholder=True),Segment('source3A','A')],mod3='OH')
    s.anneal('bottom',[Segment('samebodycomplement','xxxx',placeholder=True)],to='top',pair=('samebodycomplement','unknownbody'))
    s.overhang('top','source3A',polarity='3-prime',owner='substrate',text='SourceA-tailedDNA localendonly; fullindexedadapter/wholefragmentunassembled');return s

def protocol_graph():
    b=B('ChOR-seq published EdU/qChOR representativeH3K4me3T0 preparation',AUTHORITY,'SourceHeLaS3released2h15minafter2mM17hthymidineblock; upstreamsynchronization/outsideboundary');g=b.g
    b.op('edu_pulse','Source10uMEdU10min earlySpulse',S.MAGIC,(b.stock('edu','Actual10uMEdUnucleoside supply, notfresh invitroEdUTP',K.MOLECULAR,R.STRUCTURAL),b.cond('pulse_conditions','SourceearlyS2h15afterrelease, medium24uMdeoxycytidine')),
       reason='Positive cellularpulse; metabolicnucleotideactivation/incorporation/wholechromatinclassesunassembled. SelectednascentT0fixsimmediatelyafterpulse, no maturitychaseborrowed.')
    b.op('fix_t0','Sourceimmediate1percentformaldehydefixnascentT0; durationunassigned',S.CROSS_LINKING,(b.cond('formaldehyde','Source1percentformaldehyde'),))
    b.op('quench','Source.125Mglycine5minRT',S.QUENCHING,(b.cond('glycine','Source.125Mglycine'),))
    b.op('lyse','Sourceicecoldlysis20min/21Gneedle',S.EXTRACTION,(b.stock('lysis_preparation','Actualsource100mMNaCl/66mMTris8/5mMEDTA/.3percentSDS/1.6percentTriton/proteaseinhibitorstock, wholemixtureunassembled',K.MOLECULAR,R.STRUCTURAL),b.cond('lysis_conditions','20minonice/21Gneedle')))
    b.op('sonicate','SourceBioruptor20cycles30sON/30sOFFhigh',S.FRAGMENTATION,(b.cond('sonication','Actualsourcephysicalshearing'),))
    clear=g.product('hela_chromatin','SourceclarifiedHeLachromatinsupernatant',K.MOLECULAR);pel=g.product('hela_pellet','Actualclarificationpellet, laterfateunassigned',K.MOLECULAR)
    g.react('clarify','Source14000rpm4C10min, retainsupernatant',S.SPLIT,((b.x,R.TEMPLATE),b.cond('clarify_conditions','Sourcerpmretained, nogconversionwithoutrotor')),(clear,pel));g.mark_terminal(Terminal(pel.id,Disposition.UNRESOLVED,'Nonretainedinsolublefraction.'))
    s2=g.input('s2_chromatin','ActualasynchronousDrosophilaS2chromatinlabeled10uMEdU39h,separatelyfixed/lysed/sonicated; wholefactoryunassembled',K.MOLECULAR,primary=True,chemical_state={'source_role':'physicalexogenousreferencespike','species':'DrosophilaS2','EdU':'10uM39h'})
    mix=g.product('mixed_chromatin','ActualHeLa/S2mixedchromatin, S2.025–2.5percenttotal; sample-specificratio unassigned',K.MOLECULAR)
    g.react('spike_mix','SourcephysicalDrosophilaspike-inAFTERsonicationandBEFOREhistoneChIP',S.POOL,((clear,R.TEMPLATE),(s2,R.TEMPLATE)),(mix,));b.x=mix
    b.split('chip_input','Source30ugtotalmixedchromatinselectedforChIP','Actual30ugChIPinput','Otherpreparedchromatinstock, amount/fateunassigned')
    sn=g.product('preclear_sn','ActualProteinApreclearsupernatantkeptforantibodycapture',K.MOLECULAR);pc=g.product('preclear_bound','ActualpreclearProteinAfractionremoved; specificityunassigned',K.MOLECULAR)
    g.react('preclear','SourceProteinAagarose1h4Cpreclear',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('preclear_beads','SourceProteinAagarosepreclearpreparation',K.SELECTION,R.SELECTION),b.cond('preclear_conditions','Source4C1h;500uldialysisbuffer+400ulincubationbufferwithsourceproteaseinhibitors, exactfinalvolumebasisunassembled')),(pc,sn));g.mark_terminal(Terminal(pc.id,Disposition.UNRESOLVED,'Preclearboundfractionexcluded, notassumednonspecificpureclass.'));b.x=sn
    b.op('antibody','SourceH3K4me3C42D8antibodyovernight4C',S.BINDING,(b.stock('h3k4me3_ab','Actualsourcehistone-antibodymolecularpreparation, doseunassigned',K.MOLECULAR,R.STRUCTURAL),b.cond('antibody_conditions','Overnight4C')))
    bound=g.product('chip_bound','SourceProteinA/antibody-boundchromatinpreparation',K.MOLECULAR);flow=g.product('chip_unbound','ActualhistoneChIPnonretainedfraction',K.MOLECULAR)
    g.react('chip_capture','SourcepreblockedProteinAagarose3hcapture',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('protein_a_blocked','ActualsourceProteinAagarosepreblocked1mgmlBSAinRIPAovernight, sameblockedbeadstocknotfreshBSAaddition',K.SELECTION,R.SELECTION,chemical_state={'co_present_carrier':'Source1mgmlBSAovernightpreblocking; notjustacondition'}),b.cond('capture_conditions','Source3hincubation')),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'UnboundChIPmaterialnotsequencedinselectedroute.'));b.x=bound
    b.clean('chip_wash','SourceRIPA3x/RIPA.5MNaClonce/LiClonce/TEtwice')
    b.op('rnase','SourceRNaseA37C30min',S.DIGESTION,(b.stock('rnase_a','ActualsourceRNaseA, doseunassigned',K.CATALYST,R.CATALYST),b.cond('rnase_conditions','37C30min')))
    b.op('proteinase','Source.5percentSDS/ProteinaseK1mgml37C10h',S.DEPROTEINIZATION,(b.stock('proteinase_k','Actualsource1mgmlProteinaseK',K.CATALYST,R.CATALYST),b.cond('pk_conditions','Source.5percentSDS37C10h')))
    b.op('reverse_crosslink','Source65C6hdecrosslink',S.THERMAL_TREATMENT,(b.cond('reverse_conditions','Source65C6h'),))
    b.clean('minelute','SourceMinElutePCRDNAPurification')
    b.split('library_input','Source100ngimmunoprecipitatedDNAselectedforindexedlibrary','Actual100ngpurifiedChIPDNA','UntransferredChIPDNAstock, mass/fateunassigned')
    b.op('end_a','SourceKAPAHyperprepend-repair/A-tail beforeindexedadapterligation',S.MAGIC,(b.stock('kapa_end_preparation','ActualhistoricalKAPAend/A-tailpreparedmolecularstocks, enzymes/NTformulationsunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Sourceend-repair/A-tailpositive; local3Aclassretained, fullDNAenddistributions/historicalstockcompositionunassembled.')
    adapters=g.input('indexed_adapters','Source1.25uMPentabaseIllumina-compatibleindexedadapterstock; wholeordering/indexlength/valueunassembled',K.ADAPTER,chemical_state={'sample_index_feature':feature('chor.sample','sample_index','unknown'),'whole_adapter':'unassembled'})
    geom=AdapterLigationGeometry('unavailable','unavailable',(source_a_end(),unknown_scene()),'Known sourceA-tailedsubstrate3A retained; fullPentabaseindexedadapterstrands/ends/wholejunctionunassembled, notabsence.')
    b.op('indexed_ligate','Sourceindexedadapterligation20C60min',S.LIGATION,((adapters,R.STRUCTURAL),b.stock('kapa_ligase','ActualsourceKAPAligationcatalystcomponent, doseunassigned',K.CATALYST,R.CATALYST),b.cond('ligation_conditions','Source20C60min/1.25uMadapter')),geometry=geom)
    b.clean('postligation','SourceAMPureXPcleanupbeforeqChORpool')
    own=b.split('own_pool_input','ActualselectedT0indexedDNAcontributionforqChORtimecoursepool','SelectedT0indexedDNAaliquot','OtherT0indexedDNAstock, amount/fateunassigned')
    others=g.input('other_indexed_timepoints','ActualotherindexedDNAfromT1/T6/T12sameH3K4me3timecourse; eachseparateupstreampreparation, amounts/molecularstatesunassembled',K.MOLECULAR)
    pool=g.product('indexed_timecourse_pool','ActualsourceqChORindexedtimecourseDNApool BEFOREclick',K.MOLECULAR);g.react('pool_timepoints','SourcephysicalindexedDNApoolacrosssametimecoursebeforeEdUclick',S.POOL,((own,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    b.split('click_input','Source200ngindexedmixedDNAselectedforclick','Actual200ngpooledindexedDNA','Otherindexedtimecoursepoolstock, mass/fateunassigned')
    b.op('click','SourceClickITRT30minonALREADYindexedmixedDNA',S.CHEMICAL_CONVERSION,(b.stock('biotin_teg_azide','Actual.5mMbiotinTEGazidemolecularreactant',K.MOLECULAR,R.STRUCTURAL),b.stock('click_copper','Actual.1mMCuSO4catalyst',K.CATALYST,R.CATALYST),b.cond('click_conditions','Source1xClickITbuffer/.5mMTHPTA/5mMaminoguanidine/10mMsodiumascorbate')),
       product='SourcebiotinTEG-taggedEdUandunreactedDNA/adapter/free-labelpopulations; nativechemicalEdUnotconvertedtoordinaryT')
    b.clean('postclick_size','SourceAMPureXP200–700bpselection/resuspendTE')
    captured=g.product('t1_bound','SourceT1-retainedindexedEdUDNApreparation',K.MOLECULAR);unbound=g.product('t1_unbound','ActualT1nonretainedfraction',K.MOLECULAR)
    g.react('t1_capture','SourceMyOneT1equilibratedB&W/equalvolume2xB&W,RT30minrotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('t1_beads','ActualsourceMyOneStreptavidinT1beads',K.SELECTION,R.SELECTION),b.cond('bw_conditions','Source1xB&W5mMTris7.5/.5mMEDTA/1MNaCl/.05percentTween; equalvolume2xstock')),(captured,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'NonretainedDNAisnotcopiedintoselectedlibrarybyassumption.'));b.x=captured
    b.clean('t1_wash','Source1xB&W4x/TE.05percentTween2x/10mMTrisonce; resuspendactualbeadsinwater')
    bulk=g.input('pcr_stock','ActualsourceKAPAHyperprep10–15cyclepreparedPCRstock; fullprimers/NTs/catalystidentitydosesunassembled',K.MOLECULAR)
    e=g.product('pcr_enzyme','SamepreparedPCRenzymecomponent',K.CATALYST);nt=g.product('pcr_nt','SamepreparedPCRnucleotidecomponent',K.MOLECULAR);f=g.product('pcr_forward','Samepreparedsourceforwardprimercomponent, wholeorderingunassembled',K.OLIGO);r=g.product('pcr_reverse','Samepreparedsourcereverseprimercomponent, wholeorderingunassembled',K.OLIGO);buf=g.product('pcr_buffer','SamepreparedPCRbuffer',K.CONDITION)
    g.react('pcr_components','Representactualpreparedstockcomponents, nofreshduplicateadditions',S.GENERIC,((bulk,R.TEMPLATE),),(e,nt,f,r,buf),changes_molecule=False)
    b.op('onbead_pcr','SourceKAPA10–15cyclesusingACTUALstreptavidinbeadsuspensionastemplate',S.PCR,((e,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='Wholeindexedtemplate/PCRprimersites/products/cyclechoiceunassembled. AmplifiedDNAcopy distinctfromoriginalcovalentbiotinEdUtemplateclasses.'))
    sup=g.product('pcr_supernatant','ActualPCRsupernatant retainedaftermagneticbeadseparation; carryover/yieldunassigned',K.MOLECULAR);beads=g.product('pcr_beads','Actualmagneticallyretainedbeadresidue/originaltemplateclasses; laterfateunassigned',K.MOLECULAR)
    g.react('postpcr_magnet','SourcepostPCRmagneticseparation keepsSUPERNATANTforcleanup',S.SPLIT,((b.x,R.TEMPLATE),),(sup,beads));g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'Retainedbeadresiduenotimplicitlysequenced; originaltemplatefatesnotquantified.'));b.x=sup
    b.clean('final_ampure','SourceAMPureXPcleanuponactualPCRsupernatant')
    b.split('load_aliquot','Actualsequencingloadallocation, amountunassigned','SelectedqChORloadingaliquot','Othercleanedlibrarypoolstock, amount/fateunassigned')
    lib=g.product('library','SourceqChORH3K4me3indexedtimecoursepool',K.LIBRARY);g.react('load','SourceHiSeq4000SE50orNextSeq500SE75',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceHiSeq4000SE50/NextSeq500SE75; per-samplemachine/indexmask/wholeadapter/primerlandingunassembled. S2referencespeciescountsusedcomputationallyafteritsactualphysicalspike-in. PCRduplicateremovaldoesnotimplyUMIs. Directbiotin-dUTPchromatinpull-downbeforedecrosslinkingandotherhistone/timebranchesnotmerged.'))
    g.add_source_item(SourceItem('selected_route','PublishedEdU/qChORChIPDNA indexing→timecoursepool→click/T1→onbeadPCR→PCRsupernatant; directbiotin-dUTPchromatincapturehasdifferentorder.','SourceS1p2timepointsandS2p3sizeQC; linkedChIP/click/vendororiginalsordinaryunread. No fullsource/reaction/absenceaudit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('NascentT0fixedimmediatelyafterpulse',g.steps['fix_t0'].inputs[0].material,'edu_pulse_product')
    c('SpikeactualphysicalDNAcontribution',len(g.steps['spike_mix'].inputs),2)
    c('ChIPDNApurifiedbeforeindexedlibrary',g.steps['library_input'].inputs[0].material,'minelute_kept')
    c('Timecoursepoolphysicallybeforeclick',list(g.steps).index('pool_timepoints')<list(g.steps).index('click'),True)
    c('Clickusesactual200ngpooledaliquot',g.steps['click'].inputs[0].material,'click_input_selected')
    c('PCRusesbeadboundwashedDNA',g.steps['onbead_pcr'].inputs[0].material,'t1_wash_kept')
    c('FinalcleanupusesactualPCRsupernatant',g.steps['final_ampure'].inputs[0].material,'pcr_supernatant')
    c('No directbiotindUTPbranchimport',all('biotin_dutp'not in x for x in g.materials),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2018.08.010_PMC6202308.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2018.08.010')):
        t=p.read_text();c('Publishedpoolbeforeclicksequence','indexed DNA from all time points in the same time course were then mixed together. Then, Click-IT' in t,True);c('PublishedpostPCRsupernatant','the supernatant was cleaned-up with Agencourt AMPure XP beads' in t,True)
    c('Privatepartialremainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
