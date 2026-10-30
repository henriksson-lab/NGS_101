"""Published2026 BALI, selected single-digit RNA demonstration; private incomplete graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature
from illumina import P7,TRUSEQ_READ1,TRUSEQ_READ2,NEBNEXT_UNIVERSAL_PRIMER
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
TITLE='BALI published2026 single-digit spatial RNA route'
AUTHORITY='PublishedPNAS2026PMC13214022 main129–141 andSIAppendix321–389/TableS2originalp32. SIjoint450–555 read onlyforpositivealternativeboundary; no2024preprintauthority. Twooriginals, sourceinvestigationclosed.'
READOUT='PublishedmainR1=100/R2=200,300cyclereagents; sourceRNAbranchNovaSeq. R1mapsRNA, R2decodesspatialbarcode. Completeprimerlandings/indexmask remainunassembled. PhysicalROI tags differfromsequencingsampleindexes andUMIs; sourcebarcodealignment/filtering iscomputational.'
RT_SEQUENCE='TAA'+'AATAACCCTAG'+'N'*8+'ACTG'+'T'*30+'VN'
INDEX_ARM=TRUSEQ_READ2[13:]

def inventory(name):
    sc=Scene()
    if name=='BL621':
        sc.strand(name,[Segment('upstreamphotocleavablecap','TAA',note='Source /iSpPC/ immediatelyafterthese3bases; capnotretainedafterphotocleavage'),Segment('root','AATAACCCTAG'),Segment('RTvariablewindow','N'*8,placeholder=True,feature=feature('bali_rt_n8','inline_barcode','unknown'),note='Source8N; functionalidentityunassignedinselectedscope,notassumedtheterminalindexUMI'),Segment('RTspacer','ACTG'),Segment('polyT','T'*30),Segment('anchor','VN',placeholder=True)]);return sc
    if name=='BL617':sc.strand(name,[Segment('sourceTSOhandle','T'+TRUSEQ_READ1),Segment('sourceRNAguanosines','GGG',note='rGrGrG')]);return sc
    if name in {'BL599','BL601'}:
        sc.strand(name,[Segment('overhang','GGTTATT'),Segment('leftflank','CG'),Segment('spatialcode','ATCTAGGT' if name=='BL599' else 'GGACTCGT',feature=feature('bali_roi','spatial_barcode','fixed')),Segment('rightflank','GC')]);return sc
    if name in {'BL623','BL624'}:
        sc.strand(name,[Segment('PCRarm',INDEX_ARM),Segment('terminalUMI','N'*6,placeholder=True,feature=feature('bali_terminal_umi','umi','random')),Segment('leftflank','GC'),Segment('spatialcode','ACCTAGAT' if name=='BL623' else 'ACGAGTCC',feature=feature('bali_roi','spatial_barcode','fixed')),Segment('rightflank','CG')],mod5='deSBioTEG');return sc
    if name=='universal':sc.strand(name,[Segment('sourceuniversal',NEBNEXT_UNIVERSAL_PRIMER)]);return sc
    if name in {'A012','A0XX'}:
        sc.strand(name,[Segment('P7',P7),Segment('sampleindex','TACAAG' if name=='A012' else 'N'*6,placeholder=name=='A0XX',feature=feature('bali_sample_index','sample_index','fixed' if name=='A012' else 'unknown')),Segment('sourcetruncatedR2body',TRUSEQ_READ2[:-1])]);return sc
    return None

def chemistry(name):
    if name=='BL621':return {'polymer':'DNA','sequence':RT_SEQUENCE,'photocleavable_linker_after':3,'written':'TAA/iSpPC/AATAACCCTAGNNNNNNNNACTG'+'T'*30+'VN','after_uncaging':'Source exposed5primephosphate; exactretainedproductnotconstructed'}
    if name=='BL617':return {'polymer':'DNA/RNA','rna_positions':tuple(range(len('T'+TRUSEQ_READ1),len('T'+TRUSEQ_READ1)+3))}
    if name in {'BL623','BL624'}:return {'five_prime':'/5deSBioTEG/','source':'PublishedSITableS2p32; desthiobiotinTEGnotordinarybiotin'}
    return None

def paired_index(top_name,bottom_name):
    top=inventory(top_name).strands[top_name];bottom=inventory(bottom_name).strands[bottom_name]
    sc=Scene()
    sc.strand('top',[Segment('topoverhang',top.segs[0].top),Segment('topleftflank',top.segs[1].top),Segment('topcode',top.segs[2].top,feature=top.segs[2].feature),Segment('toprightflank',top.segs[3].top)])
    sc.anneal('bottom',[Segment('bottomPCRarm',bottom.segs[0].top),Segment('bottomUMI',bottom.segs[1].top,placeholder=True,feature=bottom.segs[1].feature),Segment('bottomleftflank',bottom.segs[2].top),Segment('bottomcode',bottom.segs[3].top,feature=bottom.segs[3].feature),Segment('bottomrightflank',bottom.segs[4].top)],to='top',pair=('bottomcode','topcode'),mod5=bottom.mod5)
    sc.ligation_end('top','topoverhang',kind='cohesive',polarity='5-prime',owner='adapter')
    return sc

def unknown():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('input','SourcefreshfrozenE16.5embryonicbraincoronalcryosection; selectedRNAbranch',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,stock,product,oligos=(),retained_inputs=(),adapter_geometry=None,catalyst=None,magic_reason=None):
        nonlocal x
        enzyme=kind in {S.REVERSE_TRANSCRIPTION,S.DEPROTEINIZATION,S.PCR,S.LIGATION,S.DIGESTION}
        if catalyst is not None:enzyme=catalyst
        c=g.input(key+'_stock',stock,K.CATALYST if enzyme else K.CONDITION);ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        ins.extend(retained_inputs)
        for i,name in enumerate(oligos):ins.append((g.input(key+'_oligo_'+str(i),name,K.OLIGO,state=inventory(name),chemical_state=chemistry(name)),R.STRUCTURAL))
        kw={}
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown(),reason='Sourceprimersknown; wholeRNA/cDNA/ligationproductandlandinggeometryunassembled.')
        if kind==S.LIGATION:kw['geometry']=adapter_geometry
        if kind==S.MAGIC:kw['reason']=magic_reason or 'Source-specificcombinedoperation;actualmoleculartransitionunassembled.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        out=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(out,),**kw);x=out
    def wash(key,label):
        nonlocal x
        c=g.input(key+'_buffer',label,K.CONDITION);y=g.product(key+'_kept','Sourcephysicallyrecoveredpreparation; recovery/carryoverunquantified',K.MOLECULAR);w=g.product(key+'_waste','Actualsourcewashfractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physicalworkup,notperfectselectivity.');x=y
    def allocation(key,label,disp):
        nonlocal x
        y=g.product(key+'_production','Sourceproductionstockafteractualallocation',K.MOLECULAR);z=g.product(key+'_allocated',label,K.QC if disp==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));
        if disp is not None:g.mark_terminal(Terminal(z.id,disp,'Sourceanalyticalmaterialnotreturnedtoproduction.' if disp==Disposition.QC else 'Untransferredsourcepreparationretainedoutsidechosenreaction.'))
        x=y
    op('equilibrate','RT5min→37C5min;papboundarydry5min',S.THERMAL_TREATMENT,'Sourceambient/heat','Sourceequilibratedsection')
    wash('oct','Source500ulPBSremovesOCT')
    op('fix','500ul.5percentPFAexactly5minRT',S.CROSS_LINKING,'SourcePFA','Sourcefixedsection')
    op('quench','1.25Mglycinequickwash+5minRT',S.QUENCHING,'Sourceglycine','Sourcequenchedsection');wash('pbs','TwoPBS5minwashes')
    op('permeabilize','.5percentTritonX100/PBS15minRT',S.PERMEABILIZATION,'SourceTriton','Sourcepermeabilizedsection');wash('perm_wash','ThreePBS5minwashes')
    op('acid','.1NHCl5min',S.CHEMICAL_CONVERSION,'SourceHCl','Sourceacidtreatedsection');wash('acid_wash','PBS2quick+one5minwash')
    op('heat','PrewarmedPBS65C5min→icecoldPBS10min',S.THERMAL_TREATMENT,'Sourceheating/cooling','SourcesectionreadyforlightprotectedRT')
    op('rt_switch','SourceBL621/dNTPfirst14ulmix65C5min→ice>2minbeforeenzymeaddition;assembled50ulRT+TS42C~16hlightprotected',S.MAGIC,'Source2.5ulSSIV/10ul5xbuffer/1ulRNAsin/5ul5Mbetaine/.5ulETSSB/2.5ul100mMDTT/2.5uldNTP10mMeach/water; sourceprimersseparateinputs','SourceRNA/cDNAandcappedroot-containingmixture; sourceRT/TScombinedtransitionunassembled',('BL621','BL617'),catalyst=True,magic_reason='SourceSSIVRT/templateswitchingknown;explicitwholecopyjunctionunassembled.')
    wash('rt_wash','PBSonce/2xSSCtwice5min');op('stain','Draq5 1:1000in2xSSC10min',S.GENERIC,'SourceDraq5','Sourceimagedsection; histologyguidesROIselection');wash('stain_wash','Two2xSSC5minwashes')
    for region,pair in [('svz',('BL599','BL623')),('cortex',('BL601','BL624'))]:
        top=g.input(region+'_index_top',pair[0],K.OLIGO,state=inventory(pair[0]),chemical_state=chemistry(pair[0]));bottom=g.input(region+'_index_bottom',pair[1],K.OLIGO,state=inventory(pair[1]),chemical_state=chemistry(pair[1]))
        cond=g.input(region+'_anneal_condition','Sourceindexpair45uM2xSSC95C5min→65C5min→4C;.1Csecramps',K.CONDITION)
        ds=g.product(region+'_index_stock','Sourcepreannealedindexstock;wholewrittenstrands/12bppairedstem/7nt5primeoverhangknown',K.ADAPTER,state=paired_index(*pair))
        g.react(region+'_index_anneal','Sourceindexstockannealingbeforeligation',S.HYBRIDIZATION,((top,R.TEMPLATE),(bottom,R.STRUCTURAL),(cond,R.CONDITION)),(ds,))
        used=g.product(region+'_index_used','Sourceindexstockaliquotforfinal3uM50ulligation;exactpreparedstockvolumeunassigned',K.ADAPTER)
        rest=g.product(region+'_index_rest','Untransferredsourcepairedindexstock;exactvolume/fateunassigned',K.ADAPTER)
        g.react(region+'_index_aliquot','Sourcepairedindexstockallocatedtochosenligation',S.SPLIT,((ds,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Indexstocknotconsumedinfullwithoutsourceallocation.'))
        geometry=AdapterLigationGeometry.known('cohesive',paired_index(*pair),statement='PublishedSITableS2p32/sourceRNApairs:12bpCG+8ntcode+GCstem,top5primeGGTTATT7ntoverhang,bottom5deSBioTEG/6NUMI/PCRtail. Mainallows6or7ntinterfaces; completeRNA/cDNAsubstrateandsealedjunctionremainunassembled.')
        op(region+'_uncage','SP5confocal30mW405nmmaxpower5minselectedROI;200ul2xSSC/coverslip',S.CHEMICAL_CONVERSION,'Sourceillumination','SourceROI5primephosphateexposed; nonilluminatedmaterialnotphysicallyremoved')
        wash(region+'_uncage_wash','2xSSCquick+5minwash then1xligasebuffer5min')
        op(region+'_ligate','50ul3uMannealedindex/T4ligase100Uul25C30min',S.LIGATION,'SourceT4ligase/50mMTrisHCl7.4/10mMMgCl2/10mMDTT/1mMATP; pairannealed45uM2xSSC95C5min→65C5min→4C ramps.1Csec','SourceROIindex/UMI/handle-bearingcDNAplusunreactedclasses',retained_inputs=((used,R.STRUCTURAL),),adapter_geometry=geometry)
        wash(region+'_lig_wash','2xSSC+Draq510min→2xSSC+10percentEC5min→2xSSC5min')
    op('lysis','200ullysisbuffer60C~16h',S.MAGIC,'Sourceprinted.8Uml/RIPAbuffer; activecomponentidentity/typeunassignedinRNAparagraph','Sourcelysatecontainingbarcodedandunbarcodedmaterial;exactdissolutionchemistryunassembled',catalyst=False,magic_reason='PublishedRNAlysisparagraphgives.8Uml/RIPAwithoutenzymeidentity;sourceknownlysisboundaryopaque,noproteaseborrowedfromjointbranch.')
    wash('minelute','SourceMinEluteReactionCleanup/elute25ulEB')
    allocation('qpcr_input','Source2ulqPCRaliquot;Q5/EvaGreen/universal+A012;40cycles98C15s68C20s72C20safter98C3min;thresholdone-thirdplateau',None)
    # qPCR has its own source primer inventory; analytical rather than production output.
    qf=g.input('qpcr_forward','universal',K.OLIGO,state=inventory('universal'));qr=g.input('qpcr_reverse','A012',K.OLIGO,state=inventory('A012'));qc=g.input('qpcr_mix','Source25ul2xQ5/.25ul50uMuniversal/.5ul50uMA012/2.5ul20xEvaGreen;50ul',K.CATALYST);qout=g.product('qpcr_readout','AnalyticalqPCRsignal;notreturnedtoproduction',K.QC)
    g.react('qpcr','Sourcealiquotdeterminescycles',S.PCR,((g.materials['qpcr_input_allocated'],R.TEMPLATE),(qf,R.STRUCTURAL),(qr,R.STRUCTURAL),(qc,R.CATALYST)),(qout,),geometry=PCRGeometry.unavailable(unknown(),reason='Analyticalampliconunassembled.'));g.mark_terminal(Terminal(qout.id,Disposition.QC,'Analyticalmaterial.'))
    allocation('unused_lysate','Source13ulremainingafter2ulQCand10ulproduction; actualstoredfateunassigned',Disposition.UNRESOLVED)
    op('pcr','Source10ulproductionintoQ5/EvaGreenPCR;cyclesfromqPCRthreshold',S.PCR,'25ul2xQ5/.25ul50uMuniversal/.5ul50uMA0XX/2.5ul20xEvaGreen;listedreactionsetas50ul','SourceamplifiedspatialRNApreparation; exactampliconunassembled',('universal','A0XX'))
    wash('ampure','.8xAMPure/elute20ulEB');allocation('tapestation','SourceTapestationD5000andKAPAlibraryquantificationanalyticalallocations',Disposition.QC)
    others=g.input('other_libraries','Sourceotherindependentlypreparedlibraries;actualpoolmembers/volumesunassigned',K.LIBRARY,primary=True);pool=g.product('pooled','Physicalsequencingpool',K.MOLECULAR);g.react('pool','SourcepoolafterlibraryQC',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));lib=g.product('library','SourceBALIspatialRNAlibrary;wholegeometryunassembled',K.LIBRARY);g.react('interface','NGSpreparationidentity',S.GENERIC,((pool,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=READOUT));g.add_source_item(SourceItem('scope',AUTHORITY,'Noabsenceclaim/completechemistry/publicpage.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Onlysource2RNAregionsligated',[k for k in g.steps if k.endswith('_ligate')],['svz_ligate','cortex_ligate'])
    c('CortexuncagingfollowsSVZworkup',g.steps['cortex_uncage'].inputs[0].material,'svz_lig_wash_kept')
    c('PCRusesproductionstockafterQC',g.steps['pcr'].inputs[0].material,'unused_lysate_production')
    c('AnalyticalqPCRseparate',g.terminals['qpcr_readout'].disposition,Disposition.QC)
    c('Productionindexunassigned',inventory('A0XX').strands['A0XX'].segs[1].top,'N'*6)
    c('WholeTSOsourceTprefixandRNA3Gs',''.join(s.top for s in inventory('BL617').strands['BL617'].segs),'TACACTCTTTCCCTACACGACGCTCTTCCGATCTGGG')
    c('WholeA012sourcebodyistruncated',''.join(s.top for s in inventory('A012').strands['A012'].segs),'CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC')
    c('WholePCrootordering',RT_SEQUENCE,'TAAAATAACCCTAGNNNNNNNNACTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN')
    cut=chemistry('BL621')['photocleavable_linker_after'];c('WholePCmodifiedrootfixture',RT_SEQUENCE[:cut]+'/iSpPC/'+RT_SEQUENCE[cut:],'TAA/iSpPC/AATAACCCTAGNNNNNNNNACTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN')
    tso=''.join(seg.top for seg in inventory('BL617').strands['BL617'].segs);rnas=set(chemistry('BL617')['rna_positions']);c('WholeRNAmodifiedTSOfixture',''.join(('r'+b)if i in rnas else b for i,b in enumerate(tso)),'TACACTCTTTCCCTACACGACGCTCTTCCGATCTrGrGrG')
    for name,raw in [('BL623','/5deSBioTEG/AGACGTGTGCTCTTCCGATCTNNNNNNGCACCTAGATCG'),('BL624','/5deSBioTEG/AGACGTGTGCTCTTCCGATCTNNNNNNGCACGAGTCCCG')]:
        strand=inventory(name).strands[name];c('Wholemodifiedsourceindex '+name,'/5'+strand.mod5+'/' +''.join(seg.top for seg in strand.segs),raw)
    for region in ('svz','cortex'):
        geo=g.steps[region+'_ligate'].geometry;c('Sourceknowncohesiveindexend '+region,(geo.availability,geo.kind),('known','cohesive'))
        c('Ligationconsumespreparedindexaliquot '+region,[(i.material,i.role)for i in g.steps[region+'_ligate'].inputs if i.role==R.STRUCTURAL],[(region+'_index_used',R.STRUCTURAL)])
    c('RNAlysisdoesnotinheritjointprotease',(g.steps['lysis'].kind,g.materials['lysis_stock'].kind),(S.MAGIC,K.CONDITION))
    c('KnownSSIVretainscatalystrole',g.materials['rt_switch_stock'].kind,K.CATALYST)
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
