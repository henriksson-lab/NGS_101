"""Published PCL-seq selected fresh-frozen tissue preparation; source-qualified partial molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature
from illumina import P5,P7,TRUSEQ_READ2
from nextera import ME,ME_RC
from rt import SMART_HANDLE
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
TITLE='PCL-seq selected fresh-frozen tissue route'
AUTHORITY='PublishedPMC12401048 mainFig2/Methods159–201; publishedSITable2originalp10 verified. Twooriginalscappedpass; RNA/cDNA/RT/S5P7geometryunassembled. FFPEpretreatmentnotmixedintoFFroute.'
READOUT='SourceNovaSeqatNovogene; actualread/indexmasksunassigned. RT7NsourceUMI retaineddistinctfromadjacentV/terminalVN anchor and7ntadapterROI/samplecodes. Sequencefilters/cellbarcodecorrectionnotphysicalDNAremoval.'
RT_PREFIX='AGGCCAGAGCATTCG'
LINKER='CGAAUGCUCUGGCCUCTCAAGCACGTGGAT'
BLOCK='ATCCACGTGCTTGAG'
ADAPTER_ARM=SMART_HANDLE+TRUSEQ_READ2

def seq(name):
    return {'rt':RT_PREFIX+'N'*7+'V'+'T'*24+'VN','linker':LINKER,'block':BLOCK,'adapter1':ADAPTER_ARM+'AAACCGG'+BLOCK,'adapter2':ADAPTER_ARM+'CAAGCCT'+BLOCK,'tso':SMART_HANDLE+'GAAT'+'GGG','ispcr':SMART_HANDLE,'me_long':'TCGTCGGCAGCGTCAG'+ME,'me_short':ME_RC,'p7':P7+'N'*6+TRUSEQ_READ2,'s5':P5+'N'*8+'TCGTCGGCAGCGTC'}.get(name)

def chemistry(name):
    if name=='rt':return {'polymer':'DNA','five_prime':'PC linker','scope':'SourcePCcapidentity;exactchemicalorderingunassembled','source_photoproduct':'MainFig2 exposes5primephosphate','umi_start':len(RT_PREFIX),'umi_length':7}
    if name in {'adapter1','adapter2'}:return {'polymer':'DNA','five_prime':'Biotin'}
    if name=='linker':return {'written':LINKER,'polymer':'Source mixedU/T ordering retained; detailed sugar partition unassigned','explicit_u_positions':tuple(i for i,b in enumerate(LINKER)if b=='U')}
    if name=='tso':return {'polymer':'DNA/RNA','rna_positions':tuple(range(len(SMART_HANDLE)+4,len(SMART_HANDLE)+7))}
    if name=='me_short':return {'polymer':'DNA','five_prime':'P'}
    return None

def inventory(name):
    text=seq(name)
    if not text:return None
    sc=Scene();parts=[]
    if name=='rt':parts=[Segment('sourceRThandle',RT_PREFIX),Segment('source7NUMI','N'*7,placeholder=True,feature=feature('pcl_umi','umi','random')),Segment('sourceprepolyTV','V',placeholder=True,note='Meaningunassigned; notsilentlyaddedto7NsourceUMI'),Segment('polyT','T'*24),Segment('anchor','VN',placeholder=True)]
    elif name.startswith('adapter'):parts=[Segment('sourcehandles',ADAPTER_ARM),Segment('sourceROI_sample_code','AAACCGG'if name=='adapter1'else'CAAGCCT',feature=feature('pcl_roi_sample','sample_index','fixed'),note='Table2distinct7ntbarcodesto distinguish samples/replicates; notanIlluminaI1index'),Segment('linkerarm',BLOCK)]
    elif name in {'p7','s5'}:
        n=len(P7)if name=='p7'else len(P5);w=6 if name=='p7'else 8
        parts=[Segment('flowcellarm',text[:n]),Segment('sampleindex',text[n:n+w],placeholder=True,feature=feature('pcl_'+name+'_index','sample_index','unknown')),Segment('sourceprimerbody',text[n+w:])]
    elif name=='tso':parts=[Segment('sourceTSOhandle',SMART_HANDLE+'GAAT'),Segment('sourceRNA3Gs','GGG',note='rGrGrG')]
    else:parts=[Segment(name,text,note='Sourceorderingonly; exactsugar/end andassembledgeometrynotinferred')]
    mods=chemistry(name)or{};sc.strand(name,parts,mod5=mods.get('five_prime',''));return sc

SOURCE_ORDERINGS={'adapter_ligate_oligo_1':{'source':'PublishedSupplementalTable2 originalPDFp10, redAAACCGGsamplebarcode','written':'/Biotin/AAGCAGTGGTATCAACGCAGAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAACCGGATCCACGTGCTTGAG'}}
SOURCE_OLIGO_EXCEPTIONS={'adapter_ligate_oligo_1':'PublishedTable2originalp10 wholebiotinadapter1verified; extractioninterleavesrownamebetweenhandleandredbarcode/linkercontinuation. Localwholeorderingfixturepreserveseveryfixedbase/biotin;notanabsenceclaim.'}

def local_adapter_interface():
    """Source Table2 local splint/adapter/cDNA-handle interface, not full product."""
    source=inventory('adapter1').strands['adapter1']
    sc=Scene()
    sc.strand('linker',[Segment('linker_rt_arm',LINKER[:15],note='Literal sourceU/T retained; RNA-equivalent complement to RT_PREFIX. Full sugar partition unassigned.'),Segment('linker_adapter_arm',LINKER[15:])],label='Source Linker')
    sc.anneal('cdna',[Segment('cdna_handle',RT_PREFIX),Segment('cdna_remainder','[remaining cDNA unassembled]',placeholder=True)],to='linker',pair=('cdna_handle','linker_rt_arm'),mod5='P',label='Primer-derived cDNA 5prime handle')
    sc.anneal('adapter',[Segment('adapter_handles',source.segs[0].top),Segment('adapter_sample_code',source.segs[1].top,feature=source.segs[1].feature),Segment('adapter_block',source.segs[2].top)],to='linker',pair=('adapter_block','linker_adapter_arm'),mod5=source.mod5,label='Source Biotin adapter1')
    # Coordinates come from the two independently supplied source strands.
    # Adjacent handle and adapter termini are a local nick, not a sealed product.
    if sc.strands['cdna'].span('cdna_handle')[1] != sc.strands['adapter'].span('adapter_block')[0]:
        raise ValueError('Source splint arms must align the cDNA5prime and adapter3prime termini at one nick')
    sc.ligation_end('adapter','adapter_block',kind='splinted',polarity='3-prime',owner='adapter',text='Source adapter3prime terminus beside cDNA5prime phosphate; product unassembled')
    sc.ligation_end('cdna','cdna_handle',kind='splinted',polarity='5-prime',owner='substrate',text='Source illuminated primer-derived cDNA5prime phosphate')
    return sc

def unknown():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[unassembled]',placeholder=True)]);return sc

def tn5_inventory():
    sc=Scene();sc.strand('loaded_long',[Segment('handle','TCGTCGGCAGCGTCAG'),Segment('mosaic',ME)])
    sc.anneal('loaded_short',[Segment('mosaic_rc',ME_RC)],to='loaded_long',pair=('mosaic_rc','mosaic'),mod5='P')
    return sc

# Literal curator-reviewed selected published preparation, not graph enumeration.
REVIEWED_INPUTS = '''input fix_stock quench_stock permeabilize_stock acid_stock neutralize_stock rt_primer rt_hyb_buffer rt_wash_buffer rt_stock block_stock block_oligo_0 block_oligo_1 block_wash_buffer uncage_stock pc_wash_buffer adapter_ligate_stock adapter_ligate_oligo_0 adapter_ligate_oligo_1 formamide_buffer initial_release_stock slide_rinse_buffer rinse_collection overnight_stock first_ampure_buffer beads capture_buffer capture_wash_buffer second_rt_ts_stock second_rt_ts_oligo_0 wta_stock wta_oligo_0 wta_clean_buffer dilute_stock loaded_tn5a_stock stop_stock tag_clean_buffer index_pcr_stock index_pcr_oligo_0 index_pcr_oligo_1 index_clean_buffer'''.split()
REVIEWED_OPERATIONS = '''fix quench permeabilize acid neutralize rt_hybridize rt_primer_wash rt_components rt block block_wash uncage pc_wash adapter_ligate formamide initial_release recover_slide_washes collect_washes overnight first_ampure capture capture_wash second_rt_ts_components second_rt_ts wta_components wta wta_clean qubit dilute tn5a_components tag stop tag_clean index_pcr_components index_pcr index_clean interface library'''.split()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('input','Sourcefreshfrozentissuecryosection10um; selectedFFroute, actualRNAsequencesunassigned',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,stock,product,oligos=(),retained_inputs=()):
        nonlocal x
        enzyme=kind in {S.REVERSE_TRANSCRIPTION,S.TEMPLATE_SWITCHING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION}
        supplied_components=key in {'rt','second_rt_ts','wta','index_pcr'}
        c=tag_enzyme if key=='tag' else g.input(key+'_stock',stock,K.MOLECULAR if supplied_components else K.CATALYST if enzyme else K.CONDITION);ins=[(x,R.TEMPLATE)]
        if supplied_components:
            cat=g.product(key+'_enzyme','Maxima H Minus component of source supplies'if key in {'rt','second_rt_ts'}else'KAPA HiFi polymerase component of source amplification supplies',K.CATALYST)
            nts=g.product(key+'_nts','Source100uM dNTP substrates'if key=='rt'else'Source1mM dNTP substrates'if key=='second_rt_ts'else'Source printed0.3M dNTP substrates; concentration convention unresolved, no mM repair',K.MOLECULAR)
            buffer=g.product(key+'_buffer','Actual source reaction-buffer component',K.CONDITION)
            parts=[cat,nts,buffer];ins.extend(((cat,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)))
            if key in {'rt','second_rt_ts'}:
                inhibitor=g.product(key+'_inhibitor','Actual source1percent RNase inhibitor contribution',K.MOLECULAR)
                parts.append(inhibitor);ins.append((inhibitor,R.STRUCTURAL))
            g.react(key+'_components','Same source reaction supplies; no extra doses, stock-count or aggregate-packaging assumption',S.GENERIC,((c,R.TEMPLATE),),parts,changes_molecule=False)
        else:ins.append((c,R.CATALYST if enzyme else R.CONDITION))
        ins.extend(retained_inputs)
        for i,name in enumerate(oligos):ins.append((g.input(key+'_oligo_'+str(i),name,K.OLIGO,state=inventory(name),chemical_state=chemistry(name)),R.STRUCTURAL))
        kw={}
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown(),reason='Sourceprimerinventoryknown; wholebiotin-cDNA/TSO/tagmentation-productandprimerlandingunassembled.',symmetric_species=oligos==('ispcr',))
        if kind==S.LIGATION:
            kw['geometry']=(AdapterLigationGeometry.known('splinted',local_adapter_interface(),statement='Publishedmain69/Fig2 andoriginalTable2p10:literal30ntU/TLinkerjoins15ntRT_PREFIX-complementand15ntBLOCK-complementarms;whole5primeBiotinadapter1andprimer-derivedcDNA5primephosphateterminialignatanadjacentnick. FullcDNA,sealedproduct,linkersugarpartitionandbiochemicalyieldunassigned.') if key=='adapter_ligate' else AdapterLigationGeometry('unavailable','unavailable',(unknown(),),'Backgroundblockingoperationknown;fullbackgroundspecies/junctiongeometryunassembled.'))
        if kind==S.GENERIC:kw['changes_molecule']=False
        out=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(out,),**kw);x=out
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_buffer',label,K.CONDITION);y=g.product(key+'_kept','Sourcephysicallyretainedpreparation; recovery/carryoverunquantified',K.MOLECULAR);w=g.product(key+'_waste','Actualsourcewash/supernatantfraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Actualphysicalworkupnotperfectmolecularselectivity.');x=y
    op('fix','FF4percentformaldehyde10minRT',S.CROSS_LINKING,'Sourceformaldehyde','Sourcefixedsection')
    op('quench','1.25Mglycine5min',S.QUENCHING,'Sourceglycine','Sourcequenchedsection')
    op('permeabilize','.5percentTriton3min',S.PERMEABILIZATION,'SourceTriton','Sourcepermeabilizedsection')
    op('acid','.1NHCl5min',S.CHEMICAL_CONVERSION,'SourceHCl','Sourceacidtreatedsection')
    op('neutralize','1MTrisHCl8.0/10minRT',S.QUENCHING,'SourceTris','Sourceneutralizedsection')
    # Primer binds before enzyme addition and is carried through actual washes, not added twice.
    rtprimer=g.input('rt_primer','rt',K.OLIGO,state=inventory('rt'),chemical_state=chemistry('rt'));hybb=g.input('rt_hyb_buffer','Source2xSSC/10percentformamide/printed50percentdextransulfate/1percentRNaseinh/30minRT',K.CONDITION);hyb=g.product('rt_hybridized','SourceprimerhybridizedtissueRNA; exactbound/unboundclassesunassembled',K.MOLECULAR);g.react('rt_hybridize','SourceRTprimerhybridization',S.HYBRIDIZATION,((x,R.TEMPLATE),(rtprimer,R.STRUCTURAL),(hybb,R.CONDITION)),(hyb,));x=hyb
    washb=g.input('rt_wash_buffer','Three2xSSCwashes/one1xRTbufferwash',K.CONDITION)
    washed=g.product('rt_primer_wash_kept','SourcewashedtissueRNAtemplatecomponent; same-tubeannealedRTprimertrackedseparately,notaphysicalfraction',K.MOLECULAR)
    retained=g.product('retained_rt_primer','SourceannealedRTprimerretainedthroughactualwashes; same-tubecomponent,notnewaddition',K.OLIGO,state=inventory('rt'),chemical_state=chemistry('rt'))
    lost=g.product('rt_wash_waste','Actualsourcewashesincludingunboundmaterial',K.WASTE)
    g.react('rt_primer_wash','Actualwashfollowedbycomponenttrackinginonetube',S.CLEANUP,((x,R.TEMPLATE),(washb,R.CONDITION)),(washed,retained,lost));g.discard(lost,'Actualwashfractions.');x=washed
    op('rt','MainRT42C90min',S.REVERSE_TRANSCRIPTION,'Source1xRTbuffer/100uMdNTP/1percentRNaseinh/MaximaHMinus25Uul','SourceRNA-cDNAhybridswithPCcappedprimerhandles; noTSOaddedatfirstRT',retained_inputs=((retained,R.STRUCTURAL),))
    op('block','10uMlinker+10uMblock15minRTfollowedT4ligation30min37C',S.LIGATION,'SourceT4DNAligase20Uul/1xT4buffer;sourcehybridization/washesprecedeenzyme','SourcebackgrounduncagedRTspeciesblockedwithnonamplificationarms; selectivitynotperfect',('linker','block'))
    clean('block_wash','65percentformamideremovesnoncovalentsingle-strandedbindingpartners')
    op('uncage','SourceROI100percentUV10sec;slidefaceup',S.CHEMICAL_CONVERSION,'SourceAndorMosaic/NikonDMDillumination','Sourceilluminated5primephosphatecDNAplusunilluminatedcapped/blockedmaterial')
    clean('pc_wash','SourcethoroughPBSwashincluding3minrinse;cleavedcapremovedphysically')
    op('adapter_ligate','Source1uMlinker/1uMadapterhybrid15minRT;2xSSC/T4washthenligase37C30min',S.LIGATION,'SourceT4DNAligase20Uul/1xT4buffer','SourcebiotinylatedROIadapter-cDNAplusincomplete/unligatedclasses',('linker','adapter1'))
    clean('formamide','Source65percentformamidewashbeforetissuelysis')
    release_buffer=g.input('initial_release_stock','Source50mMTris/1mMEDTA/.2MNaCl/ProteinaseK1mgml55C2min',K.CATALYST)
    released=g.product('initial_lysate','Actualinitiallysatetransferredfromsection;recovery/compositionunquantified',K.MOLECULAR)
    residual=g.product('residual_section','Physicalsectionremainingafterinitiallysatecollection;notanindependentDNAinput',K.BIOLOGICAL)
    g.react('initial_release','Sourceinitialproteolysisandlysatecollection55C2min',S.DEPROTEINIZATION,((x,R.TEMPLATE),(release_buffer,R.CATALYST)),(released,residual))
    rinse_buffer=g.input('slide_rinse_buffer','SourcethreeadditionalslidewashadditionswithsamereversecrosslinksolutionincludingProteinaseK1mgml;volumes/timesunassigned',K.CATALYST)
    washes=g.product('slide_washes','Actualthreewashsolutionsrecoveredfromthesameresidualsection;notnewexternalDNAstock',K.MOLECULAR)
    slide=g.product('washed_slide','Physicalslidematerialremainingafterthreewashcollections;laterfate/contentunassigned',K.BIOLOGICAL)
    rinse_condition=g.input('rinse_collection','Sourcethreesuccessivephysicalsliderinsesandcollectionintosametube;noadditionalincubationprogrammeassigned',K.CONDITION)
    g.react('recover_slide_washes','Actualthreewashcollectionsfromretainedsource-sectionmaterial',S.EXTRACTION,((residual,R.TEMPLATE),(rinse_buffer,R.CATALYST),(rinse_condition,R.CONDITION)),(washes,slide));g.mark_terminal(Terminal(slide.id,Disposition.UNRESOLVED,'Unrecoveredslidematerialnotdeclaredemptyordiscarded.'))
    pooled=g.product('collected','Physicalinitiallysateplus3collectedslidewashsolutions',K.MOLECULAR);g.react('collect_washes','Actualsourcecollectedslidewashpool',S.POOL,((released,R.TEMPLATE),(washes,R.TEMPLATE)),(pooled,));x=pooled
    # ProteinaseK remains in the collected solution; no second external enzyme addition.
    op('overnight','SourceaddSDSfinal1percent/55Covernight',S.THERMAL_TREATMENT,'SourceSDS/heat;PKalreadycarriedincollectedmixture','Sourceovernightreversecrosslinkedlysatewithproteolysisproducts')
    clean('first_ampure','Source.6xAMPure')
    beads=g.input('beads','SourcewashedMyOneStreptavidinC1',K.SELECTION);buf=g.input('capture_buffer','Source10mMTris7.5/1mMEDTA/1MNaCl;1600rpm30minRT',K.CONDITION);bound=g.product('bound','SourcebeadassociatedbiotincDNA/RNApreparation; purityunquantified',K.MOLECULAR);unbound=g.product('unbound','Actualunretainedphysicalcapturefraction',K.WASTE);g.react('capture','Sourcebiotincapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION),(buf,R.CONDITION)),(bound,unbound));g.discard(unbound,'Physicalunretainedcapturefraction; notperfectspecificity.');x=bound
    clean('capture_wash','1xSSC/.1percentSDS15minRTonce; .1xSSC/.1percentSDS10min65Cthreetimes')
    op('second_rt_ts','SourcebeadboundcontinuedcDNAextension/templateswitching42C90min;fulljunctionunassembled',S.TEMPLATE_SWITCHING,'Source1xRTbuffer/1mMdNTP/1percentRNaseinh/MaximaHMinus2Uul','SourceTSO-extendedcDNAplusRNA/hybrid/primer/reagentcompanions;combinedRT/TSgeometryunassembled',('tso',))
    op('wta','Source23cycles98C20s67C10s72C6minafter98C3min;72C5minfinal',S.PCR,'KAPAHIFIHotStartKB2500/1xbuffer/printed.3MdNTP/.2mMMgCl2/.4uMISPCR/1UulHotStartHiFi','SourcewholecDNAampliconpreparation',('ispcr',))
    clean('wta_clean','Source.6xAMPure/elute10ulTris10mMpH8')
    q=g.product('qubit_aliquot','SourceQubitanalyticalaliquot',K.QC);ret=g.product('after_qubit','SourceproductioncDNAstockafteranalyticalallocation',K.MOLECULAR);g.react('qubit','ActualsourceQubitallocation',S.SPLIT,((x,R.TEMPLATE),),(ret,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analyticalmaterialnotreturnedtoproduction.'));x=ret
    op('dilute','SourcecDNAdiluted1ngul',S.GENERIC,'Sourcediluent/volumeunassigned','Source1ngulcDNAstock')
    loaded=g.input('loaded_tn5a_stock','SourceTn5Apreparedloadedreagent; exact protein/adaptersstoichiometryandactiveloadinggeometryunassembled',K.MOLECULAR)
    tag_enzyme=g.product('loaded_tn5a_enzyme','EnzymecomponentofthesamepreloadedTn5Astock; not a separatefreeenzymeaddition',K.CATALYST)
    tag_adapter=g.product('loaded_tn5a_adapter','DNAadaptercomponentofthesamepreloadedTn5Astock; notnewfreeoligoaddition',K.ADAPTER,state=tn5_inventory(),chemical_state={'long':seq('me_long'),'short':seq('me_short'),'short_five_prime':'P','scope':'SourceTable2orderedinventory;completeactivetransposomegeometryunassembled'})
    g.react('tn5a_components','RepresentcontributingcomponentsoftheonepreloadedsourceTn5Astock; notphysicalseparationorassembly',S.GENERIC,((loaded,R.TEMPLATE),),(tag_enzyme,tag_adapter),changes_molecule=False)
    op('tag','SourceTn5A55C10min',S.TAGMENTATION,'SourceTn5A;exactactiveloadinggeometryunassembled','SourcefragmentedNexteraRead1-bearingcDNAmix;endo-gapsealedproductsunassembled',retained_inputs=((tag_adapter,R.STRUCTURAL),))
    op('stop','Source.2percentSDSquenchesTn5',S.QUENCHING,'SourceSDS','Sourcequenchedtagmentationmixture')
    clean('tag_clean','Source1xAMPure')
    op('index_pcr','Source72C5minthen15cycles95C20s55C30s72C30s;72C5minfinal',S.PCR,'KAPAHIFIKK2602/1xbuffer/printed.3MdNTP/.2mMMgCl2/.2uMP7/.2uMS5/1UulHiFi','Sourceindexedfragmentlibrary',('p7','s5'))
    clean('index_clean','SourceAMPure,ratio/elutionunassigned')
    lib=g.product('library','SourcePCL-seqNGSlibrary;wholeprimer/readgeometryunassembled',K.LIBRARY);g.react('interface','PhysicalNGSlibraryidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=READOUT));g.add_source_item(SourceItem('scope',AUTHORITY,'Wholeoligosverified;reactiongeometryandlinked-sourceauditincomplete.Noabsenceclaim.'))
    g.cover_source_item(SourceCoverage('scope','modelled','library'))
    g.declare_reviewed_scope(REVIEWED_INPUTS, AUTHORITY+' Selected actual molecular supplies and prepared-stock boundaries; unknown whole oligos remain unassembled.')
    g.declare_reviewed_scope(REVIEWED_OPERATIONS, AUTHORITY+' Actual selected physical preparation, allocations, QC and sequencing contribution; no assertion of complete occupied products or yield.')
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('RTconsumesproducedprimercomponent',[(i.material,i.role)for i in g.steps['rt'].inputs if i.role in {R.STRUCTURAL,R.SELF_PRIMER} and g.materials[i.material].kind==K.OLIGO],[('retained_rt_primer',R.STRUCTURAL)])
    c('RTstockintroducedonlybeforehybridization',sum(i.material=='rt_primer'for st in g.steps.values()for i in st.inputs),1)
    c('BlockprecedesROIillumination',g.steps['uncage'].inputs[0].material,'block_wash_kept')
    c('Formamideremovalbeforetissuelysis',g.steps['initial_release'].inputs[0].material,'formamide_kept')
    c('OvernightincubationdoesnotaddProteinaseKagain',[(i.material,i.role)for i in g.steps['overnight'].inputs if i.role==R.CATALYST],[])
    c('RinseDNAcomesfromsameresidualsection',g.steps['recover_slide_washes'].inputs[0].material,'residual_section')
    c('Remainingwashedslidephysicalfateexplicit',g.terminals['washed_slide'].disposition,Disposition.UNRESOLVED)
    c('SecondRTfollowsactualphysicalbeadwash',g.steps['second_rt_ts'].inputs[0].material,'capture_wash_kept')
    c('SecondreactionusesTSacceptornotnewRTprimer',(g.steps['second_rt_ts'].kind,[(i.material,i.role)for i in g.steps['second_rt_ts'].inputs if i.role in {R.STRUCTURAL,R.SELF_PRIMER} and g.materials[i.material].kind==K.OLIGO]),(S.TEMPLATE_SWITCHING,[('second_rt_ts_oligo_0',R.STRUCTURAL)]))
    c('SourceLinkerUandTpreserved',LINKER,'CGAAUGCUCUGGCCUCTCAAGCACGTGGAT')
    c('WholeBiotinadapter1sourceordering','/'+inventory('adapter1').strands['adapter1'].mod5+'/'+seq('adapter1'),SOURCE_ORDERINGS['adapter_ligate_oligo_1']['written'])
    c('WholeRT-PCsourceordering','/'+inventory('rt').strands['rt'].mod5+'/'+seq('rt'),'/PC linker/AGGCCAGAGCATTCGNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTVN')
    c('WholeTSOsourceordering',seq('tso'),'AAGCAGTGGTATCAACGCAGAGTGAATGGG')
    c('WholeP7sourcepattern',seq('p7'),'CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
    c('WholeS5sourcepattern',seq('s5'),'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTC')
    text=seq('tso');rnas=set(chemistry('tso')['rna_positions']);c('WholeRNAmodifiedTSOfixture',''.join(('r'+b)if i in rnas else b for i,b in enumerate(text)),'AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG')
    short=tn5_inventory().strands['loaded_short'];c('Whole5PshortMEfixture','/PHO/'+''.join(seg.top for seg in short.segs)if short.mod5=='P'else 'wrongsourceend','/PHO/CTGTCTCTTATACACATCT')
    interface=local_adapter_interface()
    c('LocalinterfacesourceLinkerretainsliteralUandT',''.join(seg.top for seg in interface.strands['linker'].segs),'CGAAUGCUCUGGCCUCTCAAGCACGTGGAT')
    adapter=interface.strands['adapter'];c('WholelocalBiotinadapterordering','/'+adapter.mod5+'/'+''.join(seg.top for seg in adapter.segs),SOURCE_ORDERINGS['adapter_ligate_oligo_1']['written'])
    c('Localnicksourceendpolarities',[(e.owner,e.polarity,e.kind)for e in interface.ligation_ends],[('adapter','3-prime','splinted'),('substrate','5-prime','splinted')])
    c('KnownROIinterfacedoesnotclaimbackgroundgeometry',(g.steps['adapter_ligate'].geometry.availability,g.steps['block'].geometry.availability),('known','unavailable'))
    for key in ('rt','second_rt_ts','wta','index_pcr'):
        c('Actual same-supply enzyme and nucleotide functions consumed '+key,
          (g.steps[key+'_components'].inputs[0].material,[(i.material,i.role)for i in g.steps[key].inputs if i.material in {key+'_enzyme',key+'_nts'}]),
          (key+'_stock',[(key+'_enzyme',R.CATALYST),(key+'_nts',R.STRUCTURAL)]))
    c('Selected reviewed preparation passes normal publication gate',g.publication_problems(),[])
    g.input('unreviewed_addition','Regression-only unreviewed stock',K.CONDITION)
    c('New stock cannot auto-publish',any('unreviewed_addition' in p and 'source-inventory' in p for p in g.publication_problems()),True);c.report()


PAGE_FILE = 'pcl_seq.html'
def render_page():
    from page import head,info
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph().require_publishable()
    parts=[head('PCL-seq'),'<div class="wrap" data-protocol-complete="true"><h1>PCL-seq</h1>',
      '<p>Photocleavable RT primers mark RNA-derived cDNA in tissue. Illumination of a chosen region exposes a ligatable end, allowing a biotinylated adapter to select its cDNA before amplification and tagmentation.</p>',
      info('Source: <a href="https://doi.org/10.1101/gr.279906.124">Published PCL-seq study</a>, fresh-frozen Methods and Supplemental Table 2. This page selects fresh-frozen tissue; FFPE pretreatment and optional immunofluorescence experiments are separate branches.'),
      '<h2>Source concentration wording</h2>',info('Published Methods195 and199 each print 0.3 M dNTPs in the amplification recipe. The model retains that literal value without silently changing it to mM. Methods167 also prints 50% dextran sulfate in the hybridization buffer. Executed formulation and concentration conventions remain unresolved; current vendor kit formulations are not substituted.'),
      '<h2>Known primers and adapter interface</h2><p>The RT primer has seven UMI bases, a separate V position and a poly-T/VN anchor. The photocleavable cap is source-declared; its complete chemical ordering is unassembled. The biotin adapters carry seven-base sample/ROI codes. These differ from the later P7/S5 sequencing indexes.</p>']
    for name in ('rt','linker','block','adapter1','adapter2','tso','ispcr','p7','s5'):
        parts.append(panel(inventory(name).rows(),caption='Source '+name+' ordering; full occupied molecules remain unassembled.'))
    parts += [panel(local_adapter_interface().rows(),caption='Known local splinted nick: source mixed-U/T linker, biotin adapter 3′ end and illuminated cDNA 5′ phosphate. Detailed linker sugar partition and whole sealed product remain unassigned.'),
      panel(tn5_inventory().rows(),caption='Source Tn5A DNA inventory with 5′-phosphorylated short mosaic end; enzyme and DNA are components of one loaded source stock.'),
      '<h2>Preparation and reads</h2><p>Hybridize the RT primer before adding enzyme, reverse transcribe and block uncaged background handles. Illuminate the ROI, remove cleaved caps and ligate the source adapter through the linker. Collect the initial lysate and three washes from the same slide, reverse crosslinks overnight, and capture biotin-bearing cDNA. Continue extension/template switching on beads, amplify, quantify, tagment and index-amplify before cleanup and sequencing.</p>',
      info('The TSO is a template-switch acceptor, not a new RNA-priming oligo. The same preloaded Tn5 stock supplies enzyme and adapter DNA. Complete occupied RNA/cDNA products, background-blocking geometry, primer landing, index values and read masks remain unresolved. Literal U/T linker bases are preserved without assigning their full sugar partition.'),
      '<p>The publication names NovaSeq sequencing at Novogene without a complete run mask for this selected trace. Computational RNA/UMI and sample-code processing does not physically purify the library.</p>',PlanLimitationSection('library'),
      '<details><summary>Full preparation and chemistry warnings</summary>',WorkflowSection('Published fresh-frozen tissue preparation'),'</details>',
      '<p><a href="01_pcl-seq.html">Source notes and model limits</a></p></div>']
    return assemble_page(g,parts)
if __name__=='__main__':selftest()
