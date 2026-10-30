"""Source-qualified CATCH-Seq reference24plex preparation; whole genomic products partial."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
AUTHORITY='Published2014PMC4214737 mainMethods47–66/Results76–85, actualoriginalPDFp2–4; actualTableS1s005DOCX header/region6cells. Twooriginals. Standard24plex reference recipe; region6PE100 acquisitionexample, not invented exact historical sample recipes.'
T7_ORDERS=('TACTACTAATACGACTCACTATAGGGT','CCCTATAGTGAGTCGTATTAGTAGTA')

def t7_oligo(which):
    sc=Scene();sc.strand('oligo',[Segment('source T7 order',T7_ORDERS[which])]);return sc

def t7_duplex():
    a,b=T7_ORDERS
    if a[:-1]!=revcomp(b):raise ValueError('source T7 duplex core no longer pairs')
    sc=Scene();sc.strand('forward',[Segment('promoter core',a[:-1]),Segment('source3T','T')]);sc.anneal('reverse',[Segment('source reverse',b)],to='forward',pair=('source reverse','promoter core'));sc.ligation_end('forward','source3T',kind='cohesive',polarity='3-prime',owner='adapter',text='Known local3T only; source ordered5phosphate/end distributions unassigned');return sc

def unknown_scene(label='Whole source product unassembled'):
    s=Scene();s.strand('unassembled',[Segment(label,'[unassembled]',placeholder=True)]);return s

def same(x):
    return {'state':x.state,'chemical_state':x.chemical_state}

class T:
    def __init__(self,g,x):self.g=g;self.x=x
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return self.g.input(key,label,kind,**kw),role
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,ports=(),kind=S.MAGIC,geometry=None,product=None):
        y=self.g.product(key+'_product',product or 'Actual source mixture; complete occupied molecules/end distribution/recovery unassembled',K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Positive source-local operation; whole occupied molecules/historical formulations/individual products and recovery unassembled, not absence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
    def clean(self,key,label,selection=True,preserve=False,selection_profile=None):
        y=self.g.product(key+'_kept','Actual retained '+label+' fraction; recovery/selectivity unassigned',K.MOLECULAR,**(same(self.x)if preserve else{}));w=self.g.product(key+'_removed','Actual removed workup fraction',K.WASTE)
        ports=((self.x,R.TEMPLATE),self.cond(key+'_conditions',label))
        if selection:ports+=(self.stock(key+'_beads',selection_profile or 'ActualSeraMag1.8xstockallocation,1:40in2.5MNaCl/20percentPEG8000;70percentEtOHwash twice/dry10min',K.SELECTION,R.SELECTION,chemical_state={'physical_stock':'SourcebeadsinNaCl/PEGbindingbuffer; no independentextraPEGdose','selection':'Source actual retained bead-boundfraction; no perfectselectivity'}),)
        self.g.react(key,label,S.CLEANUP,ports,(y,w));self.g.discard(w,'Actual excluded sourceworkup fraction; no perfectmolecular specificity claim');self.x=y
    def take(self,key,label,rest,qc=False):
        a=self.g.product(key+'_selected',label,K.MOLECULAR,**same(self.x));r=self.g.product(key+'_rest',rest,K.MOLECULAR,**same(self.x));self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r));self.g.mark_terminal(Terminal(a.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else rest));self.x=r if qc else a;return self.x
    def components(self,key,label,items):
        st=self.g.input(key+'_inventory',label,K.MOLECULAR);xs=[self.g.product(key+'_'+name,caption,kind)for name,caption,kind,role in items];self.g.react(key+'_components','Same actual aggregate supply inventory; packaging/stockcounts unassigned unless source names mix',S.GENERIC,((st,R.TEMPLATE),),tuple(xs),changes_molecule=False);return tuple((x,item[3])for x,item in zip(xs,items))

def repair_tail(t,key):
    er=t.components(key+'_repair','ActualNEBNextEndRepairModule aggregateallocation; precisepackaging/doses/historicalformulation unassigned',(('enz','Source end-repair catalyst components',K.CATALYST,R.CATALYST),('other','Source module nucleotide/cofactor/buffer components; exactidentities/doses unassembled',K.MOLECULAR,R.STRUCTURAL)))
    t.op(key+'_end_repair','SourceNEBNextEndRepairModule',er);t.clean(key+'_repair_clean','Source1.8xSPRI afterend-repair')
    at=t.components(key+'_tail','ActualNEBNextdA-TailingModule aggregateallocation; precise stockcounts/doses unassigned',(('enz','Source dA-tailcatalyst components',K.CATALYST,R.CATALYST),('other','Source tailmodule nucleotide/cofactor/buffer components; fullformulation unassembled',K.MOLECULAR,R.STRUCTURAL)))
    t.op(key+'_a_tail','SourceNEBNextdA-TailingModule',at);t.clean(key+'_tail_clean','Source1.8xSPRI afterdA-tail')

def ligase(t,key):
    return t.components(key,'ActualNEBNextQuickLigationModule aggregateallocation; packaging/formulation/doses unassigned',(('ligase','Source ligationcatalyst class',K.CATALYST,R.CATALYST),('other','Source ligationbuffer/cofactor components; identities/doses unassembled',K.MOLECULAR,R.STRUCTURAL)))

def protocol_graph():
    g=ProtocolGraph('CATCH-Seq standard24plex reference preparation and BAC-derived RNA probes',authority=AUTHORITY,representation='published_instance')
    # Actual independent sample inputs: no physical gDNA pool is invented before separate library preparation.
    shares=[]
    for i in range(1,25):
        k='s'+str(i);dna=g.input(k+'_dna','Actual independent samplegDNA; source700ng–2.5ug range/exactsamplemass unassigned; source priorQCstock volume unassigned',K.MOLECULAR,primary=True);t=T(g,dna)
        t.take(k+'_input_qc','SourceQubitHS orPicogreen and1–3ul1percentagaroseSybrGreen analytical allocations; exactchoices/volumesunassigned','Remaining independent sampleDNA afterQC',qc=True)
        t.op(k+'_volume','Adjust selectedsampleDNA to80ulwater afterinputQC',ports=(t.stock(k+'_water','Actualwaterallocationto80ul; amountunassigned',K.MOLECULAR,R.STRUCTURAL),),kind=S.GENERIC)
        t.op(k+'_shear','CovarisE210:4C/20percentduty/intensity5/200cyclesperburst/165s',ports=(t.cond(k+'_covaris','Actual80ul source referenceafterQC; exactDNArecoveryunassigned'),),kind=S.FRAGMENTATION)
        repair_tail(t,k)
        sc=Scene();sc.strand('adapter',[Segment('source inline sample adapter','[whole nonmethylated inline adapter set unassembled]',placeholder=True,feature=feature('catch.sample.'+str(i),'sample_index','unknown',note='Source distinctinlinebarcodedadapterset per sample; width/selectedvalues/fullstrandset unassembled, not molecularUMI'))])
        inv=g.input(k+'_adapter_inventory','Actual source synthesizedoligo inventory fornonmethylatedadapterset; constituentwholeorders/speciescounts unassembled',K.OLIGO,state=sc);ad=g.product(k+'_annealed_adapter','Actual sourceannealednonmethylated inlinebarcodedadapterstock',K.ADAPTER,state=sc)
        g.react(k+'_adapter_anneal','Source oligo synthesis/annealing prepares thissample adapterset',S.MAGIC,((inv,R.TEMPLATE),),(ad,),reason='Source-specific24nonmethylatedsets, not canonicalIllumina wholeorders. Actualoligostrands/annealing/endmods unassembled.')
        before=t.x;t.x=ad;adpart=t.take(k+'_adapter_take','Actual80–100pmoladapter contribution; exactsampledoseunassigned','Untransferred preparedadapterstock');t.x=before
        geom=AdapterLigationGeometry('unavailable','unavailable',(sc,),'Sourcepositiveadapterligation, but wholecustominlineadapter ends/occupiedgenomicjunctions unassembled; no canonical shell/5P borrowed.')
        t.op(k+'_adapter_ligate','Actual selectedsampleadapterset QuickLigation',((adpart,R.STRUCTURAL),)+ligase(t,k+'_ligation'),kind=S.LIGATION,geometry=geom)
        t.clean(k+'_library_clean','Source1.8xSPRI/elute40ulwater')
        t.take(k+'_library_qc','Actual1:20dilution forQubitHS/Picogreen analytical allocation; amountunassigned','Samplelibrary stock afterQC',qc=True)
        t.take(k+'_pool_quant','ActualBioanalyzer/KAPAqPCR analytical allocation forbalancedpooling','Samplelibrary stock afterpoolQC',qc=True)
        t.take(k+'_pool_share','Selected equalmolarcontributionto24plex2ug capturepool; volumesunassigned','Untransferred independent preparedsamplestock');shares.append(t.x)
    pooled=g.product('library_pool','Actual equal24plexnonmethylatedlibrarypool,2ug referenceinput',K.MOLECULAR);g.react('sample_pool','Actual physicalpoolingAFTER24independentlibrarypreparations',S.POOL,tuple((x,R.TEMPLATE)for x in shares),(pooled,))
    # BAC source collection and its physical weighted pool remain distinct from the genomic-sample array.
    bac=g.input('bac_stocks','Actual preparedselectedCTD/RP11BACsourcecollection; clonecount/individualmaps/amountsunassembled. Source200mlculture/PureLinkHiPure preparation upstream.',K.MOLECULAR);t=T(g,bac)
    t.take('bac_identity_qc','ActualQubitHS andPCR or1ugHindIII/1.5percentgelidentityQC; exactbranch/cardinalityunassigned','ActualsourceBACstocksremainingafterQC',qc=True)
    t.take('bac_pool_shares','ActualselectedweightedBACstockcontributions; clonecardinality/individualmassallocationsunassembled','UntransferredpreparedBACsource stocks afterselectedweightedcontributions')
    t.op('bac_pool','PhysicalBACpool weightedbytargetcoverage before shearing; typical4ugtotalreference',ports=(t.cond('bac_pool_boundary','Sourcearraynumber/individualmassallocations unassembled; this interface producesactualphysicalpool, not computationalconcatenation'),))
    t.op('bac_shear','SameCovarisprogramme as genomicDNA',ports=(t.cond('bac_shear_conditions','20percentduty/intensity5/200cyclesperburst/165s/4C'),),kind=S.FRAGMENTATION);repair_tail(t,'bac')
    a=g.input('t7_first','ActualsourcefirstT7adapterorder; orderedendmods unassigned',K.OLIGO,state=t7_oligo(0));r=g.input('t7_second','ActualsourcesecondT7adapterorder; orderedendmods unassigned',K.OLIGO,state=t7_oligo(1));ad=g.product('t7_adapter_stock','Actual sourceannealedT7adapterstock:26bpcore/3T',K.ADAPTER,state=t7_duplex(),chemical_state={'whole_orders':T7_ORDERS,'ordered_endmods':'unassigned; no5P invented'})
    g.react('t7_adapter_anneal','Source twoorderedoligos annealed',S.HYBRIDIZATION,((a,R.TEMPLATE),(r,R.STRUCTURAL)),(ad,))
    before=t.x;t.x=ad;adpart=t.take('t7_adapter_take','ActualselectedT7adapterligationcontribution; doseunassigned','Untransferred T7adapterstock');t.x=before
    t.op('bac_t7_ligate','SourceT7adapterligation replaces genomic inlineadapterset',((adpart,R.STRUCTURAL),)+ligase(t,'bac_ligation'),kind=S.LIGATION,geometry=AdapterLigationGeometry.known('cohesive',t7_duplex(),statement='Knownsource26bpcore/3T only; opposing5phosphate/enddistribution andsealedBACjunctions unassembled.'))
    t.clean('bac_lig_clean','Source1.8xSPRI afterT7ligation')
    t.op('t7f_anneal','100–200pmolT7F in35ul/PCRbuffer/1ul50mMMgCl2;95C5min/50C5min/4C',ports=(t.stock('t7f','ActualsourceT7F annealingoligo; wholeorder assignment toinitialstrand notexplicitlymapped, remainsunassembled',K.OLIGO,R.STRUCTURAL,state=unknown_scene('Whole sourceT7F unassembled')),t.cond('t7f_conditions','Source35ul includes10xPCRbuffer plus1ul50mMMgCl2; exactbufferadditionvolumeunassigned')),product='Source duplexpromoter/single-strandantisenseBACtemplatepreparation; exactstrandoccupancy/orientation/genomicintervals unassembled')
    t.clean('ivt_template_clean','Sourcecleanup/elute36ulwater',selection=False)
    aliquots=[g.product('ivt_template_'+str(i),'Actual12ulDNAtemplate forIVT'+str(i),K.MOLECULAR)for i in range(1,4)];g.react('ivt_allocate','Actual36ulstock dividedintoTHREE12ultemplates',S.SPLIT,((t.x,R.TEMPLATE),),tuple(aliquots));ivts=[]
    for i,x in enumerate(aliquots,1):
        k='ivt'+str(i);t.x=x
        kit=t.components(k+'_kit','ActualMegascriptT7 aggregatekitallocation; packaging/stockcounts/dosesunassembled',(('pol','SourceT7RNApolymerasecomponent',K.CATALYST,R.CATALYST),('nt','SourceRNA-nucleotide supplycomponents; ratios/formulationunassembled',K.MOLECULAR,R.STRUCTURAL),('buf','Sourcekitbuffercomponents; historicalformulationunassembled',K.CONDITION,R.CONDITION)))
        t.op(k+'_transcribe','SourceT7IVT1.5h withbiotin11UTP',ports=kit+(t.stock(k+'_bio_utp','Actualbiotin-11-UTP RNA-nucleotide contribution; labelingfraction/positionsunassigned, not barcode',K.MOLECULAR,R.STRUCTURAL),),product='ActualbiotinylatedRNAprobe reactionmixture plusDNA/unincorporatedNTPs/offproducts; wholeRNAsequences/labeloccupancyunassembled')
        t.op(k+'_dnase','SourceDNase treatment afterIVT',ports=(t.stock(k+'_dnase_stock','ActualDNaseallocation; dose/timeunassembled'),),kind=S.DIGESTION)
        t.op(k+'_stop','Actual1ul0.5MEDTA then75C10minDNaseinactivation',ports=(t.stock(k+'_edta','Actual1ul0.5MEDTA',K.MOLECULAR,R.STRUCTURAL),t.cond(k+'_heat','Source75C10min')),kind=S.QUENCHING);ivts.append(t.x)
    probe=g.product('ivt_pool','ActualphysicalpoolofTHREEcompletedIVTreactions',K.MOLECULAR);g.react('ivt_pool_step','PoolcompletedthreeIVTs AFTERseparateDNase/quench',S.POOL,tuple((x,R.TEMPLATE)for x in ivts),(probe,));t.x=probe
    t.take('probe_gel_qc','Actual2ulpooledprobe sample/95C2min/1.5percentgel10ugmlEtBr','UntransferredpooledIVTstock after2ulQC',qc=True)
    t.clean('nucaway','ActualNucAwaygel-filtration removesfreeNTPs',selection=False)
    t.take('probe_quant','ActualQubitRNA/Bioanalyzeranalyticalallocations; dosesunassigned','Probe stock afterQC',qc=True)
    t.take('probe_share','Actualselectedprobemass guidedby2500x theoreticaltargetyield; exacttargetmassunassigned','UntransferredpreparedRNAprobestock')
    probe_share=t.x;t.x=pooled
    t.op('cot1','Actual40foldhumanCot1DNA with2ug24plexpool',ports=(t.stock('cot1_dna','ActualhumanCot1DNAblocker allocation; sourceMethods40fold reference',K.MOLECULAR,R.STRUCTURAL),))
    t.op('speedvac','Concentratepool/Cot1to30ul',kind=S.GENERIC)
    t.op('library_denature','Source95C5min then65C5min',ports=(t.cond('library_heat','Actualthermalprogramme; heterogeneous strands/endoccupancyunassembled'),))
    ready=g.product('probe_ready','ActualRNAprobe workingsolution: source20ulreference plus1–2ulSUPERase-In wording',K.MOLECULAR);g.react('probe_prepare','Preheatprobe65C2min',S.MAGIC,((probe_share,R.TEMPLATE),t.stock('superase','Actual1–2ulSUPERase-In RNase-inhibitor forworkingprobesolution; exactdoseunassigned',K.MOLECULAR,R.STRUCTURAL),t.cond('probe_working_volume','Sourceprobe brought20ulwater and1–2ulinhibitor added; finalvolumeconventionunassigned')),(ready,),reason='ActualsourceRNAprobe/inhibitorworkingstock; exactactiveRNA/labeloccupancy/formulationunassembled.')
    t.op('hybridize','Source24plex104ulreference65C24–60h',ports=((ready,R.STRUCTURAL),t.stock('hyb_buffer','Actual52ul2xhybridizationbuffer10xSSPE/10xDenhardt/10mMEDTA/.2percentSDS preheated',K.MOLECULAR,R.STRUCTURAL)),product='ActualRNA-probe/genomiclibrary/Cot1hybridizationmixture; completegenomicpairing/occupancy/offtargets unassembled')
    raw=g.input('c1_raw','ActualC1streptavidinbeadstock contribution; source35ulper52ulreactionreference, exact104ulscaling/executionunassigned',K.MOLECULAR);washed=g.product('c1_washed','Actualthree-timesprewashedC1beads inbinding/washbuffer',K.SELECTION);w=g.product('c1_prewash_removed','ActualremovedC1prewashfractions',K.WASTE);g.react('c1_prewash','ThreeprewashesBEFOREhybridizationcapture;200ulbinding/washbufferper52ulreference',S.CLEANUP,((raw,R.TEMPLATE),t.cond('c1_prewash_conditions','Actualsourcebinding/washbuffer; referencevs104ulvolumeconventionunassigned')),(washed,w));g.discard(w,'ActualexcludedC1prewashfractions')
    bound=g.product('capture_product','ActualC1-retained probe/librarycomplexes; wholeoccupancy/selectivityunassembled',K.MOLECULAR);unbound=g.product('capture_unbound','Actualexcludedunboundhybridizationfraction',K.WASTE)
    g.react('capture','SourceC1capture30min/frequentpulsevortexing',S.AFFINITY_CAPTURE,((t.x,R.TEMPLATE),(washed,R.SELECTION)),(bound,unbound));g.discard(unbound,'Actualunboundcapturefractionremovedbeforebeadwashes, not perfecttargetpurity');t.x=bound
    t.clean('capture_wash','Two0.5ml1xSSC/.1percentSDSRT washes; four0.5ml.1xSSC/.1percentSDS65C washes per52ulreference',selection=False)
    eluate=g.product('standard_eluate','ActualNaOH-elutedstandardcapturedDNAlibraryfraction',K.MOLECULAR);beads=g.product('post_elution_beads','ActualresidualC1/probe/beadfraction; composition/fateunassigned',K.MOLECULAR);g.react('standard_elute','0.1MNaOH/pulsevortex10min; standardroute',S.MAGIC,((t.x,R.TEMPLATE),t.stock('naoh','Actual0.1MNaOHelutionallocation; volumeunassigned',K.MOLECULAR,R.STRUCTURAL)),(eluate,beads),reason='Sourceactualphysicalelution; wholeRNA/probe/DNAalkalineproducts andbindingyield unassembled, not inventedbiotin-cleavage.');g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'Residual bead/probe stock outside selected DNAsequencing'));t.x=eluate
    t.op('neutralize','Transfereluateto70ul1MTris7.5',ports=(t.stock('tris','Actual70ul1MTrisHCl7.5',K.MOLECULAR,R.STRUCTURAL),t.cond('neutralize_conditions','ActualTrisneutralization; finalpH/bufferconventionunassigned')),kind=S.QUENCHING)
    t.clean('capture_clean','ActualSPRI/two0.5ml70percentEtOH/dry10min/elute35ulwater',selection_profile='ActualpostcaptureSPRIallocation; ratio/stockformulationunassigned')
    t.take('pcr_template','SelectedcapturedDNAallocationfor50ulPCR; exactvolumeunassigned','Untransferred35ulcaptureeluatestock')
    primers=g.input('pcr_primer_inventory','ActualstandardIlluminaprimer supplies aggregate; wholeorders/stockcountsunassembled',K.MOLECULAR);fw=g.product('pcr_forward','Actual0.5ul25uMforwardfunctionalprimerclass; wholeorderunassembled',K.OLIGO,state=unknown_scene());rv=g.product('pcr_reverse','Actual0.5ul25uMreversefunctionalprimerclass; wholeorderunassembled',K.OLIGO,state=unknown_scene());g.react('pcr_primer_classes','Expose sourcefunctionalprimerpair, notcanonicalwholeorders orassumedpremix',S.GENERIC,((primers,R.TEMPLATE),),(fw,rv),changes_molecule=False)
    t.op('postcapture_pcr','50ulPCR/98C1min/18–22x95C30s62C3min30s',ports=((fw,R.STRUCTURAL),(rv,R.STRUCTURAL),t.stock('platinum','Actual1ulPlatinumTaq'),t.stock('pcr_dntp','Actual2.5ul10mMdNTPmix',K.MOLECULAR,R.STRUCTURAL),t.stock('betaine','Actual5ul5MBetaine',K.MOLECULAR,R.STRUCTURAL),t.cond('pcr_conditions','Actual5ul10xPCRbuffer/2ul50mMMgCl2/waterto50ul')),kind=S.PCR,geometry=PCRGeometry.unavailable(unknown_scene(),reason='ActualfunctionalstandardIlluminaprimerpair/sourcePCR; wholecustomadapter/readlanding/template moleculesunassembled; no canonicalprimerimport.'))
    t.clean('postpcr_clean','SourceSPRIafterpostcapturePCR',selection_profile='ActualpostPCRSPRIallocation; ratio/stockformulationunassigned')
    t.take('library_qc','ActualBioanalyzerHS/KAPAqPCRanalyticalallocations','Finalcapturedlibrarypool afterQC',qc=True)
    t.take('loading','Actualselectedfinalpooledloadingaliquot; concentration/volumeunassigned','UntransferredQCedcapturedlibrarypool')
    lib=g.product('library','ActualselectedstandardCATCH-Seqloadingmixture; completeproducts/readlandingunassembled',K.LIBRARY);g.react('load','SourceTableS1region6HiSeqPE100 acquisitionexample',S.GENERIC,((t.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':100,'R2':100},limitation='PublishedTableS1region6 actual24plex/HiSeqPE100. Genericreferencechemistrydoesnotassign everyrange/dose to historicalindividualsamples; wholeinlineindexes/PCR/readprimerlandingunassembled.'))
    g.add_source_item(SourceItem('scope',AUTHORITY,'Private source-richpartialgraph; not fullmolecule/assaypurity/absence validation.'))
    g.add_source_item(SourceItem('cot1_scope','Methods40foldCot1 recipe; Results typically20:1 and2.5–20series.','SelectMethodsreference40fold, preserveResultsconditions separately, not publication-errorclaim orsilent40→20repair.'))
    g.add_source_item(SourceItem('geometry_scope','Known27/26T7adapterorders yield26bpcore/3T.','Orderedendmods/T7Ffullorderassignment/occupiedBACorDNAjunctions unassembled; no inferred5P orwholecanonicalIlluminalibrary.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph();t7_duplex()
    c('Known adapter core computed from whole orders',T7_ORDERS[0][:-1],revcomp(T7_ORDERS[1]))
    c('Known adapter terminalT typed3prime',[(e.owner,e.polarity)for e in t7_duplex().ligation_ends],[('adapter','3-prime')])
    c('No invented sourceordered5P',g.materials['t7_adapter_stock'].chemical_state['ordered_endmods'],'unassigned; no5P invented')
    c('Sample pool built from24independent preparations',len(g.steps['sample_pool'].inputs),24)
    c('Every source sample is an independent primary input',all(g.materials['s'+str(i)+'_dna'].primary for i in range(1,25)),True)
    c('All actual sample contributions selected not whole-stock consumption',[x.material for x in g.steps['sample_pool'].inputs],['s'+str(i)+'_pool_share_selected'for i in range(1,25)])
    c('WeightedBACpoolusesselectedfractions',g.steps['bac_pool'].inputs[0].material,'bac_pool_shares_selected')
    c('Actual three template aliquots total36ul',g.steps['ivt_allocate'].outputs,('ivt_template_1','ivt_template_2','ivt_template_3'))
    for i in range(1,4):
        k='ivt'+str(i);c('Distinct real IVT template '+k,g.steps[k+'_transcribe'].inputs[0].material,'ivt_template_'+str(i));c('DNase follows same IVT '+k,g.steps[k+'_dnase'].inputs[0].material,k+'_transcribe_product')
    c('IVTs pooled after actual quench',[x.material for x in g.steps['ivt_pool_step'].inputs],['ivt'+str(i)+'_stop_product'for i in range(1,4)])
    c('FreeNTPcleanup follows pooled2ulQC',g.steps['nucaway'].inputs[0].material,'probe_gel_qc_rest')
    c('Probe share leaves remaining preparedstock','probe_share_rest'in g.terminals,True)
    c('ActualC1prewash producedcapturebeads',g.steps['capture'].inputs[1].material,'c1_washed')
    c('Standardeluate neutralized not bisulfitemerged',g.steps['neutralize'].inputs[0].material,'standard_eluate')
    c('PCR uses selected capturealiquot',g.steps['postcapture_pcr'].inputs[0].material,'pcr_template_selected')
    c('Unknown standardPCR pair remains typed',g.steps['postcapture_pcr'].kind,S.PCR)
    c('FinalcaptureSPRIratio not inheritedfrominitiallibraries','ratio/stockformulationunassigned'in g.materials['capture_clean_beads'].label,True)
    c('Finalloading preserves pooled remainder',g.steps['load'].inputs[0].material,'loading_selected')
    c('SourceTableS1actualPE100',g.sequencing['library'].read_lengths,{'R1':100,'R2':100})
    root=sources_root()/Path(__file__).resolve().parents[1].name
    p=root/'journal.pone.0111756_PMC4214737.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pone.0111756')):
        t=p.read_text();c('Original complete forward order','TAC TAC TAA TAC GAC TCA CTA TAG GGT'in t,True);c('Original complete reverse order','CCC TAT AGT GAG TCG TAT TAG TAG TA'in t,True)
        for title,needle in [('ActualthreeIVTs','three in vitro transcription (IVT) reactions each containing 12 ul'),('Actual2ulpooledQC','a 2 ul aliquot was mixed with gel loading buffer'),('SourceMethodsCot1','40-fold human Cot-1 DNA'),('SeparateResultsCot1','20∶1 concentration of Cot-1')]:c(title,needle in t,True)
    p=root/'journal.pone.0111756_supp_pone.0111756.s005.docx'
    if have(c,Source(p,'https://doi.org/10.1371/journal.pone.0111756.s005')):
        with zipfile.ZipFile(p)as z:
            root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};rows=[[''.join(cell.itertext())for cell in row.findall('w:tc',ns)]for row in root.findall('.//w:tr',ns)];row=next(row for row in rows if row[0]=='region 6j');c('OriginalTableS1region6plex/run row',row[:6],['region 6j','3.62','24','24','144','HiSeq 2x100bp'])
    c.report()
if __name__=='__main__':selftest()
