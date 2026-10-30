"""Private source-owned IPODinterface versusactualinput/RNAP preparations."""
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

AUTHORITY='PublishedPMC8282354main258–310/326–328; published2022PMC9045836correction; actualS5TableXLSXheaders/rows1–13. Three originals.'
ORIGINAL_CLARIFICATION_G=169000
CORRECTED_CLARIFICATION_G=15900
PHASE_SEPARATION_G=21130

def protocol_graph():
    b=B('IPOD-HR selected matched-culture interface/input/RNAP preparation',AUTHORITY,'SourceexponentialculturetargetOD600.2; upstreamgrowth/geneticmanipulationoutsidegraph');g=b.g
    rna=g.product('rna_comparator','Actual200ulculturepreservedin1mlDNA/RNAShieldforseparateRNAcomparison',K.MOLECULAR);rest=g.product('culture_production','ActualremainingcultureforIPODDNApreparation',K.BIOLOGICAL)
    g.react('rna_aliquot','Sourceremove200ulRNAcomparisonaliquotBEFORErifampin',S.SPLIT,((b.x,R.TEMPLATE),),(rna,rest));g.mark_terminal(Terminal(rna.id,Disposition.UNRESOLVED,'SeparateRNAcomparisonlibraryoutsideDNAgraph, notreturnedtoproduction.'));b.x=rest
    b.op('rifampin','Source150ugmlrifampin10minbeforecrosslinking',S.GENERIC,(b.cond('rifampin_conditions','SourcepretreatmentimmobilizesinitiatingRNAP/permitscompletion; notthecomputationalRNAPsubtraction'),))
    b.op('crosslink','Source1percentformaldehyde/10mMNaPO4RT5minvigorousshaking',S.CROSS_LINKING,(b.cond('fixative','SourceFA/NaPO4'),))
    b.op('quench','Source.333MglycineRT5min',S.QUENCHING,(b.cond('glycine','Source.333Mglycine'),))
    b.clean('cell_wash','Source2x10mlicePBSwash/dry/snapfreeze-80C; actualpelletstocks/untransferredcultureallocationunassembled')
    b.op('lyse','Source600ullysis/pellet30C15minReadyLyse',S.MAGIC,(b.stock('readylyse','Source52.5kUmlReadyLyse',K.CATALYST,R.CATALYST),b.cond('lysis_conditions','Source10mMTris8/50mMNaCl/proteaseinhibitors; stationary10xdilutionbranchnotmerged')))
    b.op('sonicate','SourceBranson25percentpower3x10s/10sicepauses',S.FRAGMENTATION,(b.cond('sonicate_conditions','Sourceweticebath'),))
    b.op('digest','Source30minicedigestionto50–200bp, calibratedpersourcestudy',S.DIGESTION,(b.stock('rnase_pre','Source60ugRNaseA',K.CATALYST,R.CATALYST),b.stock('dnase','Source6ulDNaseIFisher89835',K.CATALYST,R.CATALYST),b.cond('digestion_conditions','Source5.4ul100mMMnCl2/4.5ul100mMCaCl2; notMNase')))
    b.op('digest_stop','Source50ul500mMEDTApH8stop',S.QUENCHING,(b.cond('edta','SourceEDTAstop'),))
    clarified=g.product('clarified','Sourceclarifieddigestedlysatestock; matched-culturepooling/allocationunassembled',K.MOLECULAR);pellet=g.product('clarification_pellet','Actualnon-retainedclarificationpellet; laterfateunassigned',K.MOLECULAR)
    g.react('clarify','PublishedCORRECTED15900xg10min4C; original169000xgretainedinevidence',S.SPLIT,((b.x,R.TEMPLATE),b.cond('clarify_conditions','Actual2022correction15,900xg, notphaseextraction21,130xg')),(clarified,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Actualinsolublefractionexcludedfromselectedclarifiedstock.'));b.x=clarified
    inp=g.product('input_50','Actual50ulcommoninputsample',K.MOLECULAR);remaining=g.product('clarified_remaining','Sourceclarifiedlysateaftercommoninputselection',K.MOLECULAR)
    g.react('input_aliquot','Source50ulinputremovedbeforefractionextraction',S.SPLIT,((b.x,R.TEMPLATE),),(inp,remaining))
    ipod=g.product('ipod_lysate','ActualIPODextractionallocationofsourcepreparedstock, volumeunassembled',K.MOLECULAR);rnap=g.product('rnap_lysate','ActualmatchedRNAPChIPallocation, volumeunassembled',K.MOLECULAR);rem=g.product('lysate_remainder','Untransferredsourcepreparedlysates; quantity/fateunassigned',K.MOLECULAR)
    g.react('matched_allocations','Sourcewhenpossiblematchedculturelysatesmixedbeforeinput; exactpooling/branchquotas unassembled',S.MAGIC,((remaining,R.TEMPLATE),),(ipod,rnap,rem),reason='Sourcepossiblematched-cultureinput/IPOD/RNAPworkflow retained. Wholepelletbatchpool/aliquotallocationsequenceunassembled; nottwocopiesofthesameDNAvessel.');g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Otherpreparedstockfatesunassembled.'))
    b.x=inp;input_fraction=b.op('input_dilute','Source50ulinputdiluted1:9inelutionbuffer',S.GENERIC,(b.cond('input_elution_buffer','Source450ulelutionbuffer50mMTris8/10mMEDTA/1percentSDS'),))
    b.x=ipod;b.op('phase_mix','Source1vollysate+1vol100mMTrisbase+2vol25:24:1phenolchloroformisoamyl,vortex/RT10min',S.GENERIC,(b.cond('phase_reagents','ActualsourceTrisbase/PCIstockvolumesrelative toselectedIPODlysate'),))
    disc=g.product('white_disc','Actualretainedwhiteaqueous-organicinterface proteinDNAenrichedfraction; molecularpurityunassigned',K.MOLECULAR);aq=g.product('aqueous_phase','Actualremovedaqueousphase',K.WASTE);org=g.product('organic_phase','Actualremovedorganicphase',K.WASTE)
    g.react('phase_separate','Source21130xg2minRTseparateswhiteinterfacefrombothliquidphases',S.SPLIT,((b.x,R.TEMPLATE),b.cond('phase_spin','Source21,130xg2minRT, distinctfromcorrectedclarificationforce')),(disc,aq,org));g.discard(aq,'Sourceaqueousphasephysicallyremoved, nottheIPODDNAproduct.');g.discard(org,'Sourceorganicphasephysicallyremoved; interface retained.');b.x=disc
    b.clean('interface_wash1','Source350ulTE+350ul100mMTrisbase+700ul24:1chloroformisoamyl,vortex/21130xg2min,keepinterface')
    b.clean('interface_wash2','Source700ulTE+700ul24:1chloroformisoamyl,vortex/21130xg2min,keepinterface/wickresidualliquid')
    interface_fraction=b.op('interface_resuspend','Sourcewhiteinterfaceresuspend500ulelutionbuffer',S.GENERIC,(b.cond('interface_elution','Source50mMTris8/10mMEDTA/1percentSDS'),))
    b.x=rnap;b.op('rnap_ip_buffer','SourcemixdigestedRNAPlysate1:1with2xIPbuffer',S.GENERIC,(b.cond('ip_buffer','Source200mMTris8/600mMNaCl/4percentTriton/proteaseinhibitors2xstock'),))
    b.op('anti_rnap','Source10ulantiEcoliRNAPNeoCloneWP023overnight4C',S.BINDING,(b.stock('anti_rnap_stock','SourceactualRNAPantibodymolecularstock',K.MOLECULAR,R.STRUCTURAL),b.cond('antibody_conditions','Overnightrocking4C')))
    bead=g.product('rnap_bound','SourceProteinG/RNAPantibodyretainedfraction; capture/yieldunassigned',K.MOLECULAR);flow=g.product('rnap_flowthrough','ActualRNAPChIPnonretainedfraction',K.MOLECULAR)
    g.react('rnap_capture','Source50ulequilibratedProteinGdynabeads2h4C',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('protein_g','Actual50ulProteinGdynabeads',K.SELECTION,R.SELECTION),b.cond('rnap_capture_conditions','Source2hrocking4C')),(bead,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'NonretainedRNAPfractionnotinIPODproduction.'));b.x=bead
    b.clean('rnap_washes','Source1mleachwashA100mMTris250mMLiCl2percentTriton/B500mMNaCl1percentTriton.1deoxycholate/C500mMNaCl1percentTriton/TE; RTresuspend/reseparate')
    rnap_fraction=b.op('rnap_elute','Source500ulelutionbuffer65C30min,vortex5–10minintervals',S.MAGIC,(b.cond('rnap_elution','Source50mMTris8/10mMEDTA/1percentSDS'),),reason='Actualsourceantigenelution; fullsolid/liquidmolecularclassesandcarryoverunassembled.')
    for name,start in [('input',input_fraction),('ipod',interface_fraction),('rnap',rnap_fraction)]:
        b.x=start
        b.op(name+'_reverse','Source65C6–16hcrosslinkreversal',S.THERMAL_TREATMENT,(b.cond(name+'_reverse_conditions','Actualsourcecommoncrosslinkreversal'),))
        b.op(name+'_rnase','Source100ugRNaseA37C2h',S.DIGESTION,(b.stock(name+'_rnase_stock','Actual100ugRNaseA',K.CATALYST,R.CATALYST),b.cond(name+'_rnase_conditions','Source37C2h')))
        b.op(name+'_pk','Source200ugProteinaseK50C2h',S.DEPROTEINIZATION,(b.stock(name+'_pk_stock','Actual200ugProteinaseK',K.CATALYST,R.CATALYST),b.cond(name+'_pk_conditions','Source50C2h')))
        b.clean(name+'_dna_recover','SourcePCI/ethanolprecipitation/NaClsalt/ice95percentethanolwash',extra=(b.stock(name+'_glycoblue','ActualGlycoBluecoprecipitantmolecularcarrier, quantityunassigned',K.MOLECULAR,R.STRUCTURAL),))
        q=g.product(name+'_qc','Actualfluorescentquantitationaliquot(PicoGreenorQuantIT),amountunassigned',K.MOLECULAR);keep=g.product(name+'_dna_stock','ActualremainingproductionDNAstock',K.MOLECULAR)
        g.react(name+'_quant','SourceanalyticalDNAquantitation; optional2percentgelallocationsunassembled',S.SPLIT,((b.x,R.TEMPLATE),),(q,keep));g.mark_terminal(Terminal(q.id,Disposition.QC,'Actualanalyticalaliquotnotreturnedtoproduction.'));b.x=keep
        b.op(name+'_library_prep','SourcehistoricalNEBNextUltraDNAE7370/defaultindexkitwithsourcecleanupsubstitutions',S.MAGIC,(b.stock(name+'_library_stock','ActualhistoricalDNAkit/singleordualindexpreparedmolecularstocks; fullends/primers/NTs/catalystsdosesunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Sourcekit anddeviationspositive:prePCRZymoCleanConcentrator5/OligoCleanConcentratorinsteadAMPure; finalAMPure/Axygen, optionalrepeat. Entireinternalrepair/A-tail/ligation/PCRsequence/allocationsunassembled, nomodernconstructorinherited.')
        b.split(name+'_load_aliquot','Sourcesequencingloadingallocation, amountunassigned','Actualloadingaliquot','UntransferredsourceindexedDNAstock, fateunassigned')
        lib=g.product(name+'_library','Source'+name+'matchedDNAlibraryreadout',K.LIBRARY);g.react(name+'_load','SourceHiSeqorNextSeqpairedreadsequencing; exactsamplemaskunassigned',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
        g.sequence(lib,SequencingPlan(limitation='S5lists76x76/39x38/44x44/102x102andkitvariants; exactper-samplemask/index/wholeprimerlandingunassembled. No physicalpoolorRNAPsubtractioninvented. Separateinput/IPOD/RNAPDNAmeasurements; RNAP-adjustedIPOD-HRtrackscomputedaftersequencing.'))
    g.add_source_item(SourceItem('published_correction','Originalmain274pre-interface169000xg replacedbypublished2022correction15900xg10min4C; modelusescorrectedforce. Phaseextraction21130xg2mindistinct.','SourceMethods326–328 RNAPsignalsubtraction/inputnormalizationcomputationalonly, notphysicalRNApolymeraseorDNAdepletion.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('RNAcomparisonremovedbeforerifampin',g.steps['rifampin'].inputs[0].material,'culture_production')
    c('Inputremovedbeforeinterfaceallocation',g.steps['matched_allocations'].inputs[0].material,'clarified_remaining')
    c('Interfacewashretainsactualdisc',g.steps['interface_wash1'].inputs[0].material,'white_disc')
    c('Aqueousnotsequencingproduct',g.terminals['aqueous_phase'].disposition,Disposition.DISCARD)
    c('RNAPindependentcaptureinput',g.steps['rnap_ip_buffer'].inputs[0].material,'rnap_lysate')
    c('ThreeactualDNAlibraryclasses',set(g.sequencing),{'input_library','ipod_library','rnap_library'})
    c('No inventedRNAPphysicalsubtraction',all('subtract'not in k for k in g.steps),True)
    c('Clarificationandphasespinnotconflated',CORRECTED_CLARIFICATION_G!=PHASE_SEPARATION_G,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'journal.pbio.3001557_PMC9045836.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pbio.3001557')):
        t=p.read_text();c('Actualpublishedoldforcevalue','169,000' in t,True);c('Actualpublishedcorrectedforcevalue','15,900' in t,True)
    c('Privatepartialremainsunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
