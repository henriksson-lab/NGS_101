"""Source-qualified CATCH-Seq reference24plex preparation; whole genomic products partial."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceCoverage,SequencingPlan
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
        if selection:ports+=(self.stock(key+'_beads',selection_profile or 'ActualSeraMag1.8xstockallocation,1:40in2.5MNaCl/20percentPEG8000;70percentEtOHwash twice/dry10min',K.SELECTION,R.SELECTION,chemical_state={'physical_stock':selection_profile or 'SourcebeadsinNaCl/PEGbindingbuffer; no independentextraPEGdose','selection':'Source actual retained bead-boundfraction; no perfectselectivity'}),)
        self.g.react(key,label,S.CLEANUP,ports,(y,w));self.g.discard(w,'Actual excluded sourceworkup fraction; no perfectmolecular specificity claim');self.x=y
    def take(self,key,label,rest,qc=False):
        a=self.g.product(key+'_selected',label,K.QC if qc else K.MOLECULAR,**same(self.x));r=self.g.product(key+'_rest',rest,K.MOLECULAR,**same(self.x));self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(a,r));self.g.mark_terminal(Terminal(a.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else rest));self.x=r if qc else a;return self.x
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
    declare_page_scope(g)
    g.add_source_item(SourceItem('scope',AUTHORITY,'Source-qualified partial graph; not fullmolecule/assaypurity/absence validation.'))
    g.add_source_item(SourceItem('cot1_scope','Methods40foldCot1 recipe; Results typically20:1 and2.5–20series.','SelectMethodsreference40fold, preserveResultsconditions separately, not publication-errorclaim orsilent40→20repair.'))
    g.add_source_item(SourceItem('geometry_scope','Known27/26T7adapterorders yield26bpcore/3T.','Orderedendmods/T7Ffullorderassignment/occupiedBACorDNAjunctions unassembled; no inferred5P orwholecanonicalIlluminalibrary.'))
    g.cover_source_item(SourceCoverage('scope','modelled','library'))
    g.cover_source_item(SourceCoverage('cot1_scope','modelled','cot1'))
    g.cover_source_item(SourceCoverage('geometry_scope','modelled','t7_adapter_anneal'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as ET
    c=Check();g=protocol_graph();t7_duplex()
    c('Typed publication gate',g.publication_problems(),[])
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
    c('Postcapture selection state does not inherit initial PEG recipe',g.materials['capture_clean_beads'].chemical_state['physical_stock'],'ActualpostcaptureSPRIallocation; ratio/stockformulationunassigned')
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

# Curator-authored inventories parameterized ONLY by the source's actual24
# independent samples and THREE IVTs. Never enumerate graph keys to cover them.
def declare_page_scope(g):
    sample_targets=(
        'dna','input_qc','water','volume','covaris','shear','repair_inventory',
        'repair_components','end_repair','repair_clean_conditions','repair_clean_beads',
        'repair_clean','tail_inventory','tail_components','a_tail','tail_clean_conditions',
        'tail_clean_beads','tail_clean','adapter_inventory','adapter_anneal','adapter_take',
        'ligation_inventory','ligation_components','adapter_ligate','library_clean_conditions',
        'library_clean_beads','library_clean','library_qc','pool_quant','pool_share',
    )
    for sample in (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24):
        g.declare_reviewed_scope(tuple(f's{sample}_{target}'for target in sample_targets),
                                'Published main47–50 individual sample library; standard24plex main60/TableS1region6; '+AUTHORITY)
    g.declare_reviewed_scope((
        'sample_pool','bac_stocks','bac_identity_qc','bac_pool_shares','bac_pool_boundary',
        'bac_pool','bac_shear_conditions','bac_shear','bac_repair_inventory',
        'bac_repair_components','bac_end_repair','bac_repair_clean_conditions',
        'bac_repair_clean_beads','bac_repair_clean','bac_tail_inventory','bac_tail_components',
        'bac_a_tail','bac_tail_clean_conditions','bac_tail_clean_beads','bac_tail_clean',
        't7_first','t7_second','t7_adapter_anneal','t7_adapter_take','bac_ligation_inventory',
        'bac_ligation_components','bac_t7_ligate','bac_lig_clean_conditions',
        'bac_lig_clean_beads','bac_lig_clean','t7f','t7f_conditions','t7f_anneal',
        'ivt_template_clean_conditions','ivt_template_clean','ivt_allocate',
        'ivt_pool_step','probe_gel_qc','nucaway_conditions','nucaway','probe_quant','probe_share',
        'cot1_dna','cot1','speedvac','library_heat','library_denature','superase',
        'probe_working_volume','probe_prepare','hyb_buffer','hybridize','c1_raw',
        'c1_prewash_conditions','c1_prewash','capture','capture_wash_conditions',
        'capture_wash','naoh','standard_elute','tris','neutralize_conditions','neutralize',
        'capture_clean_conditions','capture_clean_beads','capture_clean','pcr_template',
        'pcr_primer_inventory','pcr_primer_classes','platinum','pcr_dntp','betaine',
        'pcr_conditions','postcapture_pcr','postpcr_clean_conditions','postpcr_clean_beads',
        'postpcr_clean','library_qc','loading','load','library',
    ),AUTHORITY)
    ivt_targets=('kit_inventory','kit_components','bio_utp','transcribe',
                 'dnase_stock','dnase','edta','heat','stop')
    for reaction in (1,2,3):
        g.declare_reviewed_scope(tuple(f'ivt{reaction}_{target}'for target in ivt_targets),
                                'Published main56: each of three actual12ul-template IVTs; '+AUTHORITY)

PAGE_FILE = 'catch-seq.html'

def render_page():
    from page import head,info,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    from chemistry_validation import check_graph,Status
    g=protocol_graph();g.require_publishable()
    overview=(
        '<div class="wrap" data-protocol-complete="true"><h1>CATCH-Seq</h1>'
        '<p>Make biotinylated RNA capture probes from selected BAC clones, then enrich independently prepared genomic DNA libraries. '
        'CATCH-Seq here means clone-adapted template capture hybridization; it is separate from Cas9-assisted CATCH and CATCH-IT.</p>'
        +info('Source: <a href="https://doi.org/10.1371/journal.pone.0111756">Day et al., PLOS ONE (2014)</a>, '
              'published Methods and Table S1. The model follows the standard 24-plex recipe; Table S1 region6 supplies a source-declared '
              '24-plex HiSeq paired-end100 acquisition example. It does not assign every generic recipe range to that historical cohort.')
        +caveat('The two T7-adapter orders and their local duplex are known. Whole inline sample adapters, later T7F oligo assignment, '
                'standard PCR primers, genomic/transcript sequences and complete occupied junctions remain unassembled. '
                'The bisulfite branch is not merged into this standard preparation. Passing the workflow gate is not complete molecular validation.')
        +'<h2>Two preparations meet at capture</h2>'
        +table(('Branch','Source-defined preparation','Actual allocation or pool'),(
            ('Sample libraries','QC and shear each genomic sample; repair ends, dA-tail and ligate one inline-barcode adapter set per sample.','24 independent sample inputs are prepared separately before selected balanced contributions form a 2 µg reference capture pool.'),
            ('BAC templates','QC selected BAC stocks, take weighted contributions and physically pool before shearing and T7-adapter ligation.','Unused BAC stocks remain separate; clone number and exact individual masses are unassigned.'),
            ('RNA probes','Anneal source T7F after ligation, clean up and elute 36 µl DNA template.','Three actual 12 µl template allocations feed separate T7 IVTs with biotin-11-UTP.'),
            ('Finish probes','DNase-treat and quench each IVT, then physically pool; remove 2 µl for gel QC and clean with NucAway.','Probe QC and a selected capture-probe share leave untransferred RNA-probe stock.'),
            ('Capture','Hybridize probe and pooled libraries, bind streptavidin C1 beads and wash.','Retain the bead-associated fraction; source capture is enrichment, not perfect target purity.'),
            ('Standard release','Elute with NaOH, transfer to Tris, clean up and amplify for 18–22 cycles.','Keep the standard DNA eluate separate from residual beads and the alternative bisulfite route.'),
            ('Sequence','QC the captured library and select a final loading aliquot.','Unused library and QC fractions have explicit fates; no whole-stock loading is assumed.'),
        ))
        +'<h2>Known T7-adapter duplex</h2>'
        +'<p>The source orders are 27 and 26 nucleotides. Their computed 26-base-pair core leaves a 3′ T protrusion. '
        'Ordered 5′ phosphate/end distributions are unassigned; the model does not infer them from successful ligation. '
        'The later source-named T7F oligo remains a separately declared, whole-unassembled input.</p>'
    )
    molecules=''.join(panel(t7_oligo(i).rows(),caption='Whole source T7-adapter order '+str(i+1)+'.')for i in (0,1))
    molecules+=panel(t7_duplex().rows(),caption='Computed local T7-adapter pairing and 3′ T; whole BAC junction and opposing end chemistry remain unassembled.')
    readout=('<h2>Capture conditions and read setup</h2>'
        +table(('Item','Source-defined scope'),(
            ('Reference hybridization','24-plex standard library pool, 2 µg; nominal 104 µl reaction, 65°C for 24–60 hours.'),
            ('Cot-1 wording','Methods prescribe 40-fold Cot-1. Results describe a typical 20:1 ratio and a 2.5–20-fold comparison series. The model selects the Methods reference and preserves these different scopes.'),
            ('Bead scaling','The source gives 35 µl C1 beads per 52 µl reference reaction. Exact execution/scaling for the selected 104 µl workflow remains unassigned.'),
            ('Readout example','Table S1 region6: 24-plex capture and sequencing; HiSeq 2×100 bases.'),
            ('Barcode / index layout','Inline sample-barcode meaning is known. Whole custom adapter/PCR primer bodies, selected values and precise acquired index/barcode cycle placement remain unassembled.'),
        ))
        +caveat('The nominal hybridization recipe retains its probe-volume and inhibitor-addition wording without silently normalizing volumes. '
                'Reference doses, historical sample-specific optimizations and the differing Cot-1 Results conditions are not claimed to be identical.')
    )
    claims=check_graph(g)
    audit=('<h2>Molecular validation scope</h2>'+table(('Claim status','Computed count'),tuple((s.value,str(sum(c.status==s for c in claims)))for s in Status))
           +'<p>The expanded workflow has all 24 sample preparations and three IVTs, including every QC allocation, warning and terminal fate. '
           'The repeated source inventory is parameterized by those source-declared counts, not inferred by scanning graph nodes.</p>')
    return assemble_page(g,[head('CATCH-Seq — BAC-derived RNA probes and targeted DNA capture'),overview,molecules,readout,audit,
        '<details><summary>Complete 24-sample preparation, validation warnings and material fates</summary>',
        WorkflowSection('Complete source-qualified workflow'),PlanLimitationSection('library'),'</details>',
        '<p><a href="01_catch-seq.html">Research note and source boundaries</a></p></div>'])

if __name__=='__main__':selftest()
