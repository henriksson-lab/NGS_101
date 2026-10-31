"""PHOTON published selected HeLa stress-granule route; private incomplete model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition,SourceCoverage
from chemdraw import Scene,Segment
from nextera import ADAPTOR_S7,S7_RC
from illumina import P7
from rt import SMART_HANDLE
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
TITLE='PHOTON selected cell route'
AUTHORITY='PublishedPMC12078482 Methods192–222/Fig1; publishedSI MOESM1/Table1 originalp6. Two originals, capped selected pass, not a full linked-source or absence audit.'
RT_WRITTEN='/5Alex546N//iSpPC/G/iNPOM-dT/C/iNPOM-dT/CG/iNPOM-dT/GGGC/iNPOM-dT/CGGAGATGTGTATAAGAGACAG/iBiodT/'
POLYT_WRITTEN='TT TTT TTT TTT TTT TTT TTT TTV N'
TAIL=POLYT_WRITTEN.replace(' ','')
P7_TAIL='CCCATGATCGTCCGA'
SPLINT='CCGAGCCCACGAGACTCGGACGATCATGGG'

def seq(name):
    return {'rt':ADAPTOR_S7+'T'+TAIL,'random':ADAPTOR_S7+'T'+'N'*6,'splint':SPLINT,'p7handle':P7+P7_TAIL,'tso':SMART_HANDLE+'GAATGGG','pcr_tso':SMART_HANDLE,'pcr_p7':P7}.get(name)

def inventory(name):
    sc=Scene()
    if name in {'rt','random'}:
        parts=[Segment('source_RT_handle',ADAPTOR_S7,note='N-POM thymidines at positions2,4,7,12; originalTable1 preserves literal modifications'),Segment('internal_biotin_dT','T',note='/iBiodT/'),Segment('RNA_priming_tail',TAIL if name=='rt' else 'N'*6,placeholder=True,note='Priming degeneracy, not a declared UMI')]
        sc.strand(name,parts,mod5='Alexa546/photocleavable linker')
    elif name=='tso':sc.strand(name,[Segment('TSO_handle',SMART_HANDLE+'GAAT'),Segment('terminal_RNA_Gs','GGG',note='rGrGrG')])
    else:sc.strand(name,[Segment(name,seq(name))])
    return sc

def chemistry(name):
    if name in {'rt','random'}:return {'five_prime':'/5Alex546N//iSpPC/','internal_npom_dT_positions_1based':(2,4,7,12),'internal_biotin_dT_position_1based':len(ADAPTOR_S7)+1,'photoproduct':'SourceFig1/Methods exposes5primephosphate anduncagesfourTs;yieldnotasserted'}
    if name=='tso':return {'terminal':'rGrGrG'}
    return None

def local_interface():
    sc=Scene();sc.strand('splint',[Segment('RT_complement',SPLINT[:15]),Segment('adapter_complement',SPLINT[15:])])
    sc.anneal('cdna',[Segment('RT_handle',ADAPTOR_S7[:15]),Segment('remaining_cDNA','[remaining primer/cDNA unassembled]',placeholder=True)],to='splint',pair=('RT_handle','RT_complement'),mod5='P')
    sc.anneal('adapter',[Segment('P7',P7),Segment('adapter_tail',P7_TAIL)],to='splint',pair=('adapter_tail','adapter_complement'))
    if sc.strands['cdna'].span('RT_handle')[1]!=sc.strands['adapter'].span('adapter_tail')[0]:raise ValueError('Source splint must place cDNA5prime and adapter3prime at adjacent nick')
    sc.ligation_end('adapter','adapter_tail',kind='splinted',polarity='3-prime',owner='adapter',text='SourceP7handle3prime beside uncagedRT5prime phosphate')
    sc.ligation_end('cdna','RT_handle',kind='splinted',polarity='5-prime',owner='substrate',text='SourcephotocleavedRT5primeP')
    return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unknown','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('cells','SourceHeLacellsafterSGinduction; RNAidentityunassigned',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,stock,product,oligos=(),extra=()):
        nonlocal x
        cat=kind in {S.REVERSE_TRANSCRIPTION,S.TEMPLATE_SWITCHING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.DIGESTION} and key!='dash'
        physical=key in {'fix','permeabilize','dilute','xt'}
        c=g.input(key+'_stock',stock,K.MOLECULAR if physical or key in {'rt','ts','pcr5'} else K.CATALYST if cat else K.CONDITION)
        ins=[(x,R.TEMPLATE),*extra]
        if key in {'rt','ts','pcr5'}:
            enz=g.product(key+'_enzyme','Actual source enzyme component of supplied reaction mix',K.CATALYST)
            nts=g.product(key+'_nts','Actual source nucleotide substrates of same reaction mix',K.MOLECULAR)
            buffer=g.product(key+'_buffer','Actual source buffer component of same reaction mix',K.CONDITION)
            components=[enz,nts,buffer]
            if key=='rt':
                inhibitor=g.product(key+'_inhibitor','Actual0.4Uul RNase inhibitor in source RTmix',K.MOLECULAR);components.append(inhibitor);ins.append((inhibitor,R.STRUCTURAL))
            g.react(key+'_components','Copresent source reaction-mix components; no extra fresh doses',S.GENERIC,((c,R.TEMPLATE),),components,changes_molecule=False)
            ins.extend(((enz,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)))
        else:
            ins.append((c,R.STRUCTURAL if physical else R.CATALYST if cat else R.CONDITION))
            if physical and kind in {S.CROSS_LINKING,S.PERMEABILIZATION}:
                env=g.input(key+'_environment','Source exposure/incubation conditions; physical reagent stock separately inventoried',K.CONDITION);ins.append((env,R.CONDITION))
        for i,name in enumerate(oligos):ins.append((g.input(key+'_oligo_'+str(i),name,K.OLIGO,state=inventory(name),chemical_state=chemistry(name)),R.STRUCTURAL))
        kw={}
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown(),reason='WholeboundcDNA/TSOampliconandlandingunassembled;sourceprimerinventoryknown.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry.known('splinted',local_interface(),statement='SourceTable1wholeP7/splint and RT5primehandle align adjacent nick; remainingcDNA and sealedproductunassembled.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Sourcekit boundary known; complete molecular inputs, products and internal transformations not assembled locally.'
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        b=g.input(key+'_buffer',label,K.CONDITION);y=g.product(key+'_kept','Physicalretainedsourcepreparation; purity/recoveryunassigned',K.MOLECULAR);w=g.product(key+'_waste','Actualworkupwash/supernatantfraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalremovalnotperfectselectivity');x=y
    def analytical_allocation(key,label):
        nonlocal x
        assay=g.product(key+'_assay',label,K.QC)
        kept=g.product(key+'_retained','Actual unassayed production stock; aliquot volume unassigned',K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(kept,assay))
        g.mark_terminal(Terminal(assay.id,Disposition.QC,'Analytical material does not return to production'))
        x=kept
    op('fix','Source.4percentPFA10minRT',S.CROSS_LINKING,'SourcePFA','Fixedcells');clean('pbs1','Source3PBSwashes')
    op('permeabilize','Source.5percentTriton/.4UulRNaseinh15minRT',S.PERMEABILIZATION,'SourceTriton/RNaseinh','Permeabilizedcells');clean('pbs2','Source3PBSwashes')
    op('rt','Source.5uMeachRT/random/Maxima10Uul30minRTthen37Covernight/dark',S.REVERSE_TRANSCRIPTION,'SourceMaxima/dNTP1mM/RTbuffer/RNaseinh','SourcephotocagedRNA-cDNAmixture; wholeproductsunassembled',('rt','random'))
    op('uncage','SourceSG405nm2.5mW100usROIdwell',S.CHEMICAL_CONVERSION,'Source405nmlight','IlluminatedRT5P/uncagedTsandunilluminatedcappedmaterial; noselectionperfectionclaim')
    op('release','SourcePK1to50/55C30mindark',S.DEPROTEINIZATION,'SourcePK/reversecrosslinkbuffer','LysedsourceRNA/cDNApreparation');clean('extract','SourceNucleoSpinGelPCRCleanUpXS')
    spl=g.input('splint','splint',K.OLIGO,state=inventory('splint'));p7=g.input('p7handle','p7handle',K.OLIGO,state=inventory('p7handle'));buf=g.input('adapter_buffer','Source10uMeach/STE50ul95C2min→20C1Cmin',K.CONDITION);stock=g.product('adapter_stock','Source50ulannealedsplint/P7stock; bindingoccupancyunassigned',K.ADAPTER)
    g.react('anneal_adapter','Independentsourceadapterstockannealing',S.HYBRIDIZATION,((spl,R.TEMPLATE),(p7,R.STRUCTURAL),(buf,R.CONDITION)),(stock,))
    selected=g.product('adapter_10ul','Source10ulselectedannealedadapterstock',K.ADAPTER);remain=g.product('adapter_40ul','Source40uluntransferredannealedstock',K.ADAPTER);g.react('adapter_aliquot','Source10of50uladapterallocation',S.SPLIT,((stock,R.TEMPLATE),),(selected,remain));g.mark_terminal(Terminal(remain.id,Disposition.STORED,'Untransferredsourceadapterstock; lateruseunassigned'))
    op('ligation','Source30minRTpreincubationthenT4ligase2.5ul/30minRT50ul',S.LIGATION,'SourceT4DNAligase/ligasebuffer/RNaseinh','SourceP7-ligatedcDNAandunligatedcompanions',extra=((selected,R.STRUCTURAL),))
    op('dilute','Sourceadd50ulwater',S.GENERIC,'Sourcewater','Source100ulligationpreparation')
    bead=g.input('beads','Source10ulC1beadswashedthen100ul2xB&T',K.SELECTION);b=g.input('capture_buffer','Source60minRTdarkrotate',K.CONDITION);y=g.product('bound','PhysicalstreptavidinretainedbiotincDNAandcompanions',K.MOLECULAR);w=g.product('capture_unbound','Actualcaptureunretainedfraction',K.WASTE);g.react('capture','Sourcecapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(bead,R.SELECTION),(b,R.CONDITION)),(y,w));g.discard(w,'Actualremovedsupernatantnotperfectmolecularpartition');x=y;clean('beadwash','Source3B&W-Tand1STEwashes')
    op('ts','SourceMaxima10Uul/TSO2.5uM/50ul30minRTthen42C90min',S.TEMPLATE_SWITCHING,'SourceMaxima/RTbuffer/dNTP1mM','SourcecontinuedcDNAextension/TSOproductplusRNAcontext; wholejunctionunassembled',('tso',));clean('tswash','Source3STEwashes')
    op('pcr5','Source55ulQ5/5cycles98C10s65C30s72C3min',S.PCR,'SourceQ5HotStartMasterMix/.45uMprimers','Source5cycleamplificationmixture',('pcr_tso','pcr_p7'))
    qc=g.product('qpcr_2_5ul','Source2.5ulsupernatantfor10ulqPCR;35cycleprogramme',K.QC);ret=g.product('production_remainder','RemainingproductionPCRstockincludingbeads; notanalyticalreaction',K.MOLECULAR);g.react('qpcr_split','ActualsourceqPCRallocation',S.SPLIT,((x,R.TEMPLATE),),(ret,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'AnalyticalqPCR1thirdsaturationguidesproductioncycles; DNA notreturned'));x=ret
    op('pcr_extra','Sourceadditionalcyclesfromanalytical1thirdsaturation;retainedPCRmix',S.THERMAL_TREATMENT,'Sourcecontinuedthermocycling/noexternalenzymeaddition','SourceproductionWTAproduct')
    clean('pcr_recover','Sourcecollectsupernatant/magneticbeadsremovedthen.8xAMPure');
    analytical_allocation('qubit_qc','Actual source Qubit dsDNA HS aliquot after WTA recovery; volume unassigned')
    op('xt','SourceNexteraXTmanufacturerboundary',S.MAGIC,'SourceNexteraXTkit FC1311024','Sourcekitfragmentlibrary; exactloadedends/indexprimersunassembled');clean('xtclean','Source.7xSPRISelect')
    q=g.product('tapestation_qc','SourceTapeStationanalyticalallocation;volumeunassigned',K.QC);ret=g.product('pre_dash','SourceproductionlibraryafterQC',K.MOLECULAR);g.react('qc','Physicalanalyticalallocation',S.SPLIT,((x,R.TEMPLATE),),(ret,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Actualanalysisnotreturned'));x=ret
    rnp=g.input('cas9_rnp','SourcepreassembledDASHCas9/pooledcrRNA:tracrRNAstock; guidesunassembled',K.CATALYST)
    op('dash','Source20ullibrary10nM+70ulRNP+10ulbuffer37C60min',S.DIGESTION,'Source10ul10xCas9NucleaseReactionBuffer;RNPinputseparate','SourceDASHcleavedtargetandremaininglibrarymixture',extra=((rnp,R.CATALYST),))
    op('dash_release','SourcePK1ul56C10min',S.DEPROTEINIZATION,'SourceproteinaseK','SourceproteinreleasedDASHmixture');clean('dashclean','Source.8xSPRISelect')
    analytical_allocation('post_dash_qc','Actual final post-DASH TapeStation concentration/size assay; volume unassigned')
    lib=g.product('library','SourcePHOTONsequencinglibrary; completeproduct/readgeometryunassembled',K.LIBRARY);g.react('identity','PhysicalfinalNGSlibrary',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceNextSeq2000orNovaSeqXPlus; exactcyclemasksunassigned. RandomRT6NisnotdeclaredUMI. DASHsubsetroute selected; mitochondriaSEQuoiaalternative notmerged.'))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Source-knownlocaladapterinterfaceonly;wholeproducts/readoutunassembled.'));return _publication_coverage(g)

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('WholeP7handle',seq('p7handle'),'CAAGCAGAAGACGGCATACGAGATCCCATGATCGTCCGA')
    c('WholeSplint',seq('splint'),'CCGAGCCCACGAGACTCGGACGATCATGGG')
    c('WholeTSO',seq('tso'),'AAGCAGTGGTATCAACGCAGAGTGAATGGG')
    c('CagedRTbodyafterliteralmodtranslation',seq('rt'),'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGT'+POLYT_WRITTEN.replace(' ',''))
    c('Sourcecagepositions',chemistry('rt')['internal_npom_dT_positions_1based'],(2,4,7,12))
    c('Adapteraliquotremainder',g.terminals['adapter_40ul'].disposition,Disposition.STORED)
    c('Onlyproducedadapterentersligation',[(i.material,i.role)for i in g.steps['ligation'].inputs if i.role==R.STRUCTURAL],[('adapter_10ul',R.STRUCTURAL)])
    c('TStypedacceptor',g.steps['ts'].kind,S.TEMPLATE_SWITCHING)
    c('AnalyticalPCRnotreturned',g.terminals['qpcr_2_5ul'].disposition,Disposition.QC)
    c('ContinuedPCRdoesnotaddstockenzyme',[(i.material,i.role)for i in g.steps['pcr_extra'].inputs if i.role==R.CATALYST],[])
    c('Knownlocaladapterend',g.steps['ligation'].geometry.availability,'known')
    c('Distinct source preXT and postDASH analytical allocations',g.steps['qubit_qc'].inputs[0].material=='pcr_recover_kept' and g.steps['post_dash_qc'].inputs[0].material=='dashclean_kept',True)
    c('RT nucleotide class derives from actual source RT mix',g.producer_of('rt_nts').id,'rt_components')
    c('Source-qualified publication gate',g.publication_problems(),[]);c.report()


PAGE_FILE='photon.html'
SOURCE_TARGETS=(('Publishedmain192–200 selectedpreparedHeLaSG/fix/permeabilize/RT/uncaging/release andSI1Table1 originalp6wholemodifiedprimers; fixedcells are sourceinputs, not invented biologicalsequence.', ('cells', 'fix_stock', 'fix_environment', 'fix', 'pbs1_buffer', 'pbs1', 'permeabilize_stock', 'permeabilize_environment', 'permeabilize', 'pbs2_buffer', 'pbs2', 'rt_stock', 'rt_components', 'rt', 'rt_oligo_0', 'rt_oligo_1', 'uncage_stock', 'uncage', 'release_stock', 'release', 'extract_buffer', 'extract')), ('Publishedmain200–210 adapterstock10of50ul/splintedligation/C1capture/continuedTS/5cycleQ5/actualqPCRaliquot/retainedPCR; exactlocalinterface fromsameSI1Table1.', ('splint', 'p7handle', 'adapter_buffer', 'anneal_adapter', 'adapter_aliquot', 'ligation_stock', 'ligation', 'dilute_stock', 'dilute', 'beads', 'capture_buffer', 'capture', 'beadwash_buffer', 'beadwash', 'ts_stock', 'ts_components', 'ts', 'ts_oligo_0', 'tswash_buffer', 'tswash', 'pcr5_stock', 'pcr5_components', 'pcr5', 'pcr5_oligo_0', 'pcr5_oligo_1', 'qpcr_split', 'pcr_extra_stock', 'pcr_extra', 'pcr_recover_buffer', 'pcr_recover', 'qubit_qc')), ('Publishedmain212–222 NexteraXTmanufacturerinterface/cleanup/TapeStation/selectedDASHRNP/PK/cleanup/acquisition. Historical kitwholeoligos andexactrun masks unassembled, not inheritedfromkitname.', ('xt_stock', 'xt', 'xtclean_buffer', 'xtclean', 'qc', 'cas9_rnp', 'dash_stock', 'dash', 'dash_release_stock', 'dash_release', 'dashclean_buffer', 'dashclean', 'post_dash_qc', 'identity')))

def _publication_coverage(g):
    g.cover_source_item(SourceCoverage('selected','modelled','library'))
    for source,targets in SOURCE_TARGETS:
        for target in targets:
            key='reviewed_'+target
            g.add_source_item(SourceItem(key,'Selected source scope: '+target,source))
            g.cover_source_item(SourceCoverage(key,'modelled',target))
    return g.require_publishable()

def render_page():
    from page import head,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection
    g=protocol_graph()
    content=[head('PHOTON'),'<main class="wrap" data-protocol-complete="true"><h1>'+'PHOTON'+'</h1><p>'+'Photocaged reverse-transcription primers allow spatial illumination to select handles for cDNA recovery. Subsequent splinted ligation and template switching produce a source-qualified sequencing preparation from the selected regions.'+'</p><p><a href="'+'https://doi.org/10.1038/s41467-025-59801-3'+'">Defining published article</a> · <a href="01_chem.html">Research note and source scope</a></p>',caveat('<strong>Source-qualified preparation.</strong> Actual source-known inputs and process order are retained. Unknown occupied molecules, yields and whole library/read geometry remain explicit; local reagent views do not prove a unique reaction product.'),'<h2>Preparation</h2><ol>'+'<li>Start with prepared HeLa stress-granule samples. Fix, permeabilize and reverse transcribe in the dark using dual-caged oligo(dT) and random primers.</li><li>Illuminate selected regions at 405 nm, release nucleic acids and add a source-prepared splint/P7 handle. Take 10 µl from the 50 µl adapter stock.</li><li>Capture biotin-containing material on C1 beads and continue cDNA extension/template switching. Amplify five cycles; use a separate 2.5 µl analytical qPCR aliquot to choose additional cycles for the retained reaction.</li><li>Quantify recovered WTA by Qubit, then apply the historical Nextera XT interface, clean up and assess the library by TapeStation. The selected stress-granule route applies DASH, proteinase K and cleanup, followed by a separate final TapeStation assay before sequencing.</li>'+'</ol>','<h2>Known molecular parts</h2>',table(('Part','Role and source limit'),(('Caged RT primers', 'Source Alexa546/photocleavable 5′ linker, four N-POM-caged thymidines and internal biotin-dT. Photoproduct occupancy and whole copied molecules remain unassembled.'), ('Splint and P7 handle', 'Whole source oligos support a computed local nick beside the uncaged RT 5′ phosphate. This does not reconstruct the complete cDNA.'), ('Random N6', 'A source priming class, not a declared UMI.'), ('Historical finishing', 'Nextera XT and DASH source interfaces; whole kit adapters, guide sequences and exact read landing remain unassembled.')))]
    content.append(panel(local_interface().rows(),caption='Source-known splinted handle interface; remaining cDNA and sealed product unassembled'))
    content.extend(['<h2>Read setup</h2><p>'+'NextSeq 2000 or NovaSeq X Plus. Exact run cycles, index members and full read masks are unassigned. The source spatial region is a physical illumination selection, not a post-sequencing computational filter.'+'</p>','<h2>Source limits</h2><p>'+'This page starts from prepared HeLa stress-granule samples; upstream induction/imaging-probe manufacture is outside this boundary. It follows published Methods192–222 and SI1 Table1. Other tissue/ROI illumination settings and the mitochondrial SEQuoia alternative are separate source choices. Full occupied RT/TSO/kit molecules and quantitative photolysis/selection remain unresolved.'+'</p>','<details><summary>Full workflow, stock fractions and unresolved chemistry</summary>',WorkflowSection('Selected published preparation'),'</details></main>'])
    return assemble_page(g,content)

if __name__=='__main__':selftest()
