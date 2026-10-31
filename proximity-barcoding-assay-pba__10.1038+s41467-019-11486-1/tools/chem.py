"""Published PBA exosome preparation; private incomplete graph, no MPX inheritance."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition,SourceCoverage
from chemdraw import Scene,Segment,feature
from pcr import PCRGeometry
from circular import CircularizationGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
TITLE='Proximity barcoding assay: selected exosome route'
AUTHORITY='PublishedPMC6710248 mainMethods163–187/Fig1; publishedSI MOESM1 Note1/Table1 originalp11. Two selected originals; no full linked-source or absence audit.'
# Source Ab1 whole ordering is authoritative; exact feature boundaries are declared separately.
AB1='AAGCGTGTGATTCCTAGTAATGAACAATGG'+'N'*8+'TGCTATTATGATGTCTCAGGT'
PADLOCK='GCTATTATGATGTCTCAGGU'+'N'*15+'AGAUCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTAGTGCTGGATGATCGTCC'
TEMPLATE='CTGAGACATCATAATAGCGGACGATCATCCAGCACT'
BLOCK='CTAGTGCTGGATGATCGTCC'+'U'*5

def inventory(name):
    sc=Scene()
    if name=='ab1':
        # Whole source ordering known; fixed eight-base proteinTag boundary remains unassigned locally.
        parts=[Segment('source_Ab1_fixed',AB1.split('N')[0]),Segment('moleculeTag','N'*8,placeholder=True,feature=feature('pba_molecule','umi','random')),Segment('common_arm',AB1.split('N')[-1])];mod5='thiol-MC6-D'
    elif name=='padlock':parts=[Segment('padlock_start',PADLOCK[:20]),Segment('complexTag','N'*15,placeholder=True,feature=feature('pba_complex','molecular_complex_barcode','random')),Segment('padlock_rest',PADLOCK[35:])];mod5='P'
    else:parts=[Segment(name,{'template':TEMPLATE,'block':BLOCK}[name],note='Five terminalUare2primeO-methylRNA'if name=='block'else'Source ligationtemplate andRCAprimer')];mod5=''
    sc.strand(name,parts,mod5=mod5);return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unknown','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('exosomes','Sourceisolatedexosomestock; source-specificisolationupstreamunassembled',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,stock,product,extra=(),opaque=False):
        nonlocal x
        cat=kind in {S.PCR,S.EXTENSION,S.DIGESTION};c=g.input(key+'_stock',stock,K.CATALYST if cat else K.CONDITION);ins=[(x,R.TEMPLATE),(c,R.CATALYST if cat else R.CONDITION),*extra];kw={}
        if kind==S.PCR:
            for j in range(2):ins.append((g.input(key+'_primer_'+str(j),'SourcePCRprimerpairmember;exactorderingassignmentunassembled',K.OLIGO),R.STRUCTURAL))
            kw['geometry']=PCRGeometry.unavailable(unknown(),reason='Mainpba-fwd/rev andfwd-index/rev-index namesnotmappedheretoTable1Universalprimer/Dindexrows; actualwholeproductsunassembled.')
        if kind==S.MAGIC:kw['reason']='Sourceoperationpositive; fullmoleculargeometry/internalcomponentassignmentunassembled.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def wash(key,label):
        nonlocal x
        b=g.input(key+'_buffer',label,K.CONDITION);y=g.product(key+'_kept','Physicalretainedsourcepreparation; composition/recoveryunassigned',K.MOLECULAR);w=g.product(key+'_waste','Actualsourcewash/cleanupunretainedmaterial',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalworkupnotperfectmolecularselection');x=y
    # Independent reagent preparation; no biological substrate defaults borrowed.
    pad=g.input('padlock','Sourcepadlock15NcomplexTag',K.OLIGO,state=inventory('padlock'),chemical_state={'five_prime':'P','literal_U_positions_1based':(20,39),'scope':'SourceUlettersretained;polymermodificationdetailunassigned'})
    t=g.input('template','Sourceligationtemplate',K.OLIGO,state=inventory('template'));lig=g.input('circle_ligase','SourceT4ligase.01Uul/100nMeach37C30min;buffer/ATP separately inventoried',K.CATALYST);circle=g.product('circle_mix','Sourcecircularizedpadlockplusligationtemplate/unligatedcompanions; closuregeometryunassembled',K.MOLECULAR)
    g.react('circle','Sourcepadlockcircularization',S.CIRCULARIZATION,((pad,R.TEMPLATE),(t,R.STRUCTURAL),(lig,R.CATALYST),(g.input('circle_atp','Source1mM ATP in actual circularization supplies',K.MOLECULAR),R.STRUCTURAL),(g.input('circle_buffer','Source1xphi29buffer10mMMgacetate/66mMKacetate/.1percentTween/1mMDTT',K.CONDITION),R.CONDITION)),(circle,),geometry=CircularizationGeometry.known('splinted single-strand ligation',statement='Publishedmain167ligates5primephosphorylatedpadlockonsuppliedligationtemplate; sourceTable1wholeorderingsknown. Fullpairedclosure/productstateunassembled;noevent/yieldclaim.'))
    exo=g.input('exonucleases','SourceExoI.2Uul/ExoIII2Uul37C30minthen85C20min',K.CATALYST);dig=g.product('circle_enriched','Sourceexonuclease-treatedcirclemixture; actualsurvivingprimers/linearfragmentsunassigned',K.MOLECULAR);g.react('exo','Sourceunligatedpadlockdepletion',S.DIGESTION,((circle,R.TEMPLATE),(exo,R.CATALYST)),(dig,))
    adj=g.input('primer_adjustment','SourceRCAprimersameorderingasligationtemplate adjusted100nM;additionalamountunassigned',K.OLIGO,state=inventory('template'));b=g.input('reanneal_condition','Source37C20min',K.CONDITION);templ=g.product('rca_template','Sourcecirclecomponentafterreannealing',K.MOLECULAR);primer=g.product('retained_rca_primer','Co-retainedsourceRCAprimercomponent;sametube,notphysicalfraction',K.OLIGO,state=inventory('template'));g.react('reanneal','SourceRCAprimeradjustment/reannealing',S.HYBRIDIZATION,((dig,R.TEMPLATE),(adj,R.STRUCTURAL),(b,R.CONDITION)),(templ,primer))
    phi=g.input('phi29','Sourcephi29.1Uul/dAGCU1mM37C10min',K.CATALYST);rca=g.product('rca_reaction','Sourceuracil-containingRCAproductspluscompanions; repeats/geometryunassembled',K.MOLECULAR);g.react('rca','SourceRCA',S.EXTENSION,((templ,R.TEMPLATE),(primer,R.STRUCTURAL),(phi,R.CATALYST),(g.input('rca_nts','Actual d(A,U,G,C)TP1mM substrates; source uracil-containing DNA RCA',K.MOLECULAR),R.STRUCTURAL)),(rca,))
    h=g.input('rca_heat','Source65C10min',K.CONDITION);rcastock=g.product('rca_stock','Sourceheat-terminatedRCAstockkeptminus20C;totalvolumeunassigned',K.MOLECULAR);g.react('rca_stop','SourceRCAheatstop',S.THERMAL_TREATMENT,((rca,R.TEMPLATE),(h,R.CONDITION)),(rcastock,));take=g.product('rca_selected','Source25ulRCAproductat1nM;stockdilution/transferamountunassigned',K.MOLECULAR);rem=g.product('rca_remainder','UntransferredRCAstock;amountunassigned',K.MOLECULAR);g.react('rca_aliquot','SourceRCAstockallocation',S.SPLIT,((rcastock,R.TEMPLATE),),(take,rem));g.mark_terminal(Terminal(rem.id,Disposition.STORED,'Sourcekeptstock;untransferredfraction'))
    probes=g.input('probe_pool','Sourcepreconjugated40probe/38proteinantibodypool20nMeach;SMCC/reduction/dialysisupstreamunassembled',K.MOLECULAR,chemical_state={'unpositioned_protein_identifier':feature('pba.protein','feature_barcode','unknown',note='Source proteinTag meaning known; whole conjugate/tag coordinate and specimen assignment unassembled')})
    op('stain','Sourceantibody-exosomebinding4Covernight',S.BINDING,'SourcePBAbufferinclsalmonDNA/BSA','SourceexosomeswithboundandfreeAOCpool',((probes,R.STRUCTURAL),))
    selected=g.product('one_ul','Source1ulstainedstock',K.MOLECULAR);left=g.product('stain_remainder','Untransferredstainedstock;amountunassigned',K.MOLECULAR);g.react('stain_aliquot','Actualsource1ulallocation',S.SPLIT,((x,R.TEMPLATE),),(selected,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferredstocklateruseunassigned'));x=selected
    op('dilute','Source1ulcomplexesdilutedin25ulPBAbuffer aswritten',S.GENERIC,'SourcePBAbuffer','Sourcedilutedcomplexstock;25ultotalvsdiluentnotresolved')
    aliquot=g.product('capture_20ul','Source20ulcaptureinput',K.MOLECULAR);left=g.product('dilution_remainder','Sourceuntransferreddilution;volumeunassigned',K.MOLECULAR);g.react('capture_aliquot','Actualsource20ulallocation',S.SPLIT,((x,R.TEMPLATE),),(aliquot,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferreddilution'));x=aliquot
    ctb=g.input('ctb_wells','SourceSTVcoatedtubewith50ngbiotinCTBovernightthen2washes;preparedcapturestock',K.SELECTION);b=g.input('capture_condition','Source15minRT',K.CONDITION);y=g.product('captured','SourceCTBboundexosomes/AOCs;occupancy/purityunassigned',K.MOLECULAR);w=g.product('uncaptured','Physicalunretainedcapturematerial',K.WASTE);g.react('capture','Sourceexosomecapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(ctb,R.SELECTION),(b,R.CONDITION)),(y,w));g.discard(w,'Sourcecapturedpreparationnotallvesicles');x=y;wash('wash1','Source2x100ulwashingbuffer')
    block=g.input('block','Source500nMblockingoligo',K.OLIGO,state=inventory('block'),chemical_state={'terminal':'mUmUmUmUmU','scope':'2primeO-methylRNA bases;notordinarydU'})
    op('bind_rca','Source25ulRCA1nM/block500nM37C15min',S.HYBRIDIZATION,'Sourcebindingtemperature/buffer','SourceRCA-boundreporterpreparation',((take,R.STRUCTURAL),(block,R.STRUCTURAL)));wash('wash2','Source2x100ulwashingbuffer')
    # AOC primers are components already present on captured exosomes, not free additions.
    ext=g.input('t4polymerase','Source1.25UT4DNApolymerase/25ulbuffer/dNTP100uMRT15min',K.CATALYST);y=g.product('extended','SourceAOCextensioncopiescomplexTag; someunextendedoligosmayundergo3primeexonucleolysis',K.MOLECULAR);g.react('extend','SourceproximityextensiononboundRCAtemplate;fullgeometryunassembled',S.MAGIC,((x,R.TEMPLATE),(ext,R.CATALYST),(g.input('extension_nts','Actual100uM d(A,T,G,C)TP in25ulT4polymerasebuffer',K.MOLECULAR),R.STRUCTURAL)),(y,),reason='SourceboundAOC/RCAcopyingpositive; exactsubstrate/retainedAOCprimercomponentandfullcopiedproductunassembled.');x=y;wash('wash3','SourcewashbeforePCR')
    # Explicit components of each single source bulk addition, not extra additions.
    def bulk(key,label,enzyme=False):
        stock=g.input(key+'_bulk',label,K.MOLECULAR)
        f=g.product(key+'_forward','Forwardprimercomponentofsamebulkstock;exactorderingunassembled',K.OLIGO)
        r=g.product(key+'_reverse','Reverseprimercomponentofsamebulkstock;exactorderingunassembled',K.OLIGO)
        b=g.product(key+'_buffer','Actual buffer and co-present dNTP/dye/DMSO molecular inventory from same bulk addition; no fresh duplicate dose or unspecified recipe inferred',K.MOLECULAR)
        outs=[f,r,b];enz=None
        if enzyme:
            enz=g.product(key+'_enzyme','SourcePhusion.02Uul+UNG.02Uulcomponentsofsamebulkaddition',K.CATALYST);outs.append(enz)
        g.react(key+'_components','Representsamebulkstockcomponents;notphysicalseparation',S.GENERIC,((stock,R.TEMPLATE),),outs,changes_molecule=False)
        return f,r,b,enz
    f,r,b,enz=bulk('first_mix','Source30ulPhusionHF/Phusion.02Uul/UNG.02Uul/dNTP/SYBR/DMSO/pba-fwd+rev500nMeach;exactoligomappingunassembled',True)
    y=g.product('pcr_product','SourceUNG-treated/amplifiedextensionproductreactionincludingretainedenzymes;abasicRCAcompanionsnotphysicallyremoved',K.MOLECULAR)
    g.react('pcr','SourcecombinedUNG37C15min→95C2min→30PCRcycles95C15s60C30s72C30s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.STRUCTURAL)),(y,),reason='SourceoneUNG/Phusionbulkadditionandorderedprogramme;exactprimer/productassignmentunassembled.');x=y
    take=g.product('index_1ul','Source1ulwholeamplificationreactionforindexPCR',K.MOLECULAR);left=g.product('pcr_remainder','Sourceuntransferredampliconreactionstock',K.MOLECULAR);g.react('index_aliquot','Actualsource1ulallocation',S.SPLIT,((x,R.TEMPLATE),),(take,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'UntransferredPCRstock'));x=take
    def retained(key):
        nonlocal x
        template=g.product(key+'_template','DNAsubstratecomponentofselectedwholePCRreaction;notaphysicalfraction',K.MOLECULAR)
        enz=g.product(key+'_retained_enzyme','Co-transferredpriorPhusion/UNGcomponent;amountandactivityunassigned,noexternaladdition',K.CATALYST)
        g.react(key+'_carryover_components','Representco-retainedreactioncomponents;notnewstockaddition',S.GENERIC,((x,R.TEMPLATE),),(template,enz),changes_molecule=False);x=template;return enz
    enz=retained('index')
    f,r,b,_=bulk('index_mix','Source10ulindexPCRmixcontaining10uMeachfwd-index/rev-index;newpolymeraseidentity/amountnotassigned')
    y=g.product('index_product','Source2cycleindexedPCRreactionincludingco-transferredenzymes;fullproductunassembled',K.MOLECULAR,chemical_state={'unpositioned_sample_index':feature('pba.sample','sample_index','unknown',note='Source sample indexing step; actual index values, width and whole primer coordinates unassembled')})
    g.react('index','Source1ulspikedin10ulindexmix/2cycles95C15s60C60s72C60s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.STRUCTURAL)),(y,),reason='SourceindexPCRpositive;co-transferredpriorcatalystamount/activityandcompleteprimer/productsunassembled,newenzymeadditionnotinferred.');x=y
    take=g.product('library_selected','SelectedindexedPCRreactionfor20folddilution;transferamountunassigned',K.MOLECULAR);left=g.product('index_remainder','UntransferredindexedPCRreactionstock;amountunassigned',K.MOLECULAR);g.react('library_aliquot','Actualsourceindexedstockallocation',S.SPLIT,((x,R.TEMPLATE),),(take,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferredindexedstock'));x=take
    enz=retained('library')
    f,r,b,_=bulk('library_mix','Sourcefreshmix1xPhusionHFbuffer/dNTP/SYBR/DMSO/library-fwd+rev500nM;buffernameisnotafreshenzymeaddition')
    y=g.product('library_pcr_product','Source20folddilution/15cyclelibraryPCRreaction;fullproducts/catalystactivityunassembled',K.MOLECULAR)
    g.react('library_pcr','Source20folddilutionintonewmixthen15cycles95C15s60C30s72C30s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.STRUCTURAL)),(y,),reason='SourcefreshmixlistsPhusionHFbufferandprimers;newpolymeraseadditionnotinferred.Co-transferredpriorenzymeshaveunassignedamount/activity.');x=y
    # Multiple sample pooling is not expanded into fabricated other samples.
    other=g.input('other_pcr_contributions','Other actual source PCR products contributing to the library pool; identities, amounts and upstream preparations unassigned',K.MOLECULAR)
    pooled=g.product('pool_product','Actual pool of selected product and other study PCR contributions; membership and proportions unassigned',K.MOLECULAR)
    g.react('pool_boundary','Source actual PCR product pooling before QIAquick; no cloning of one product into several samples',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,));x=pooled
    wash('purify','SourceQIAquickPCRpurification')
    lib=g.product('library','SourcePBAsequencinglibrary;wholeprimerlayoutunassembled',K.LIBRARY);g.react('library_identity','Physicallibraryidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMiSeqv2 300cyclesorNextSeq75SE;perrunexactread/indexmasksunassigned. Bioinformaticminimumreadfiltersarecomputational,notDNAremoval.'))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Physicalexosome/RCAtemplatecopyingknown;exactcircle/extension/readgeometryunassembled.'));return _publication_coverage(g)

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('WholepublishedAb1',AB1,'AAGCGTGTGATTCCTAGTAATGAACAATGGNNNNNNNNTGCTATTATGATGTCTCAGGT')
    c('Wholepublishedpadlock',PADLOCK,'GCTATTATGATGTCTCAGGUNNNNNNNNNNNNNNNAGAUCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTAGTGCTGGATGATCGTCC')
    c('Wholepublishedtemplate',TEMPLATE,'CTGAGACATCATAATAGCGGACGATCATCCAGCACT')
    c('WholeoriginalblockingRNA',BLOCK,'CTAGTGCTGGATGATCGTCCUUUUU')
    c('RCAconsumesproducedprimer','retained_rca_primer' in {i.material for i in g.steps['rca'].inputs if i.role in {R.STRUCTURAL,R.SELF_PRIMER}},True)
    c('RCAstockremainderexplicit',g.terminals['rca_remainder'].disposition,Disposition.STORED)
    c('Captureusesproduced20ul',g.steps['capture'].inputs[0].material,'capture_20ul')
    c('Sourcebindthenwashbeforeextension',g.steps['extend'].inputs[0].material,'wash2_kept')
    c('UNGPhusionaddedonesourcemix',[(i.material,i.role)for i in g.steps['pcr'].inputs if i.role==R.CATALYST],[('first_mix_enzyme',R.CATALYST)])
    c('Indexconsumesretainednotnewcatalyst',[(i.material,i.role)for i in g.steps['index'].inputs if i.role==R.CATALYST],[('index_retained_enzyme',R.CATALYST)])
    c('Libraryconsumesretainednotbufferasnewenzyme',[(i.material,i.role)for i in g.steps['library_pcr'].inputs if i.role==R.CATALYST],[('library_retained_enzyme',R.CATALYST)])
    c('Firstbulkprimerscontributemolecules',{i.material for i in g.steps['pcr'].inputs if i.role==R.STRUCTURAL} >= {'first_mix_forward','first_mix_reverse','first_mix_buffer'},True)
    c('Indexedstockremainderexplicit',g.terminals['index_remainder'].disposition,Disposition.STORED)
    c('Source1ulindexingstockremainder',g.terminals['pcr_remainder'].disposition,Disposition.STORED)
    c('Physicalselectionnotreadfilter',g.steps['purify'].kind,S.CLEANUP)
    c('Actual pool consumes independent source contribution',g.steps['pool_boundary'].kind==S.POOL and 'other_pcr_contributions' in {i.material for i in g.steps['pool_boundary'].inputs},True)
    c('Source RCA consumes nucleotide substrates','rca_nts' in {i.material for i in g.steps['rca'].inputs},True)
    c('Source-qualified publication gate',g.publication_problems(),[]);c.report()


PAGE_FILE='pba.html'
SOURCE_TARGETS=(('Publishedmain167/RCApreparation; sameSI1Note1/Table1originalp11padlock/template. Circleproductgeometry/primeramounts andRCArepeatoccupancyunassembled.', ('padlock', 'template', 'circle_ligase', 'circle_atp', 'circle_buffer', 'circle', 'exonucleases', 'exo', 'primer_adjustment', 'reanneal_condition', 'reanneal', 'phi29', 'rca_nts', 'rca', 'rca_heat', 'rca_stop', 'rca_aliquot')), ('Publishedmain173–179 selectedpreparedisolatedexosome/preconjugatedAOCpool/CTBcapture/RCAblock/proximityextension/washes. Upstreamisolation/AOCmanufacture preparedinputboundaries; noSTV-onlystarredrecipeinheritance.', ('exosomes', 'probe_pool', 'stain_stock', 'stain', 'stain_aliquot', 'dilute_stock', 'dilute', 'capture_aliquot', 'ctb_wells', 'capture_condition', 'capture', 'wash1_buffer', 'wash1', 'block', 'bind_rca_stock', 'bind_rca', 'wash2_buffer', 'wash2', 't4polymerase', 'extension_nts', 'extend', 'wash3_buffer', 'wash3')), ('Publishedmain183–187 actualone30ulUNG/PhusionPCRmix/1ulindexallocation/2cycles/retainedenzymeclass/20folddilution15cyclePCR/actualpool/QIAquick/acquisition. Freshenzyme additions or exactprimer-to-tableassignments notinvented.', ('first_mix_bulk', 'first_mix_components', 'pcr', 'index_aliquot', 'index_carryover_components', 'index_mix_bulk', 'index_mix_components', 'index', 'library_aliquot', 'library_carryover_components', 'library_mix_bulk', 'library_mix_components', 'library_pcr', 'other_pcr_contributions', 'pool_boundary', 'purify_buffer', 'purify', 'library_identity')))

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
    content=[head('Proximity barcoding assay'),'<main class="wrap" data-protocol-complete="true"><h1>'+'Proximity barcoding assay'+'</h1><p>'+'Antibody–DNA reporters on captured exosomes copy a shared tag from an RCA product. Sequencing these reporter tags links protein measurements to the same physical complex; it does not directly sequence an exosome genome.'+'</p><p><a href="'+'https://doi.org/10.1038/s41467-019-11486-1'+'">Defining published article</a> · <a href="01_chem.html">Research note and source scope</a></p>',caveat('<strong>Source-qualified preparation.</strong> Actual source-known inputs and process order are retained. Unknown occupied molecules, yields and whole library/read geometry remain explicit; local reagent views do not prove a unique reaction product.'),'<h2>Preparation</h2><ol>'+'<li>Prepare the source padlock/RCA reagent, retaining its primer and uracil-containing DNA components. Allocate the selected RCA contribution and keep its stock remainder separate.</li><li>Begin the assay with isolated exosomes and a prepared antibody–DNA pool. Stain, allocate the source 1 µl contribution, dilute and take the actual 20 µl capture input for CTB-coated wells.</li><li>Bind RCA products and blocking oligo, wash, and extend captured reporters with T4 DNA polymerase. The blocker limits extension into adjacent RCA monomers; complete extension products remain unassembled.</li><li>Add the source UNG/Phusion mixture. Allocate 1 µl of amplification reaction for indexing, then dilute indexed product into the library PCR. Carryover enzyme amount/activity remains unassigned; no extra fresh enzyme is invented.</li><li>Pool the selected PCR product with other actual source PCR contributions and purify with QIAquick before sequencing.</li>'+'</ol>','<h2>Known molecular parts</h2>',table(('Part','Role and source limit'),(('Reporter tags', 'Source eight-base molecule tags are molecular identifiers. The fixed protein-tag position and actual pool members remain unassigned in the whole conjugate.'), ('RCA complex tag', 'Fifteen random bases identify a source complex. The RCA template is physically copied onto reporters, not a computationally assigned pixel.'), ('Blocking oligo', 'Whole source sequence with five terminal 2′-O-methyl U residues. These are modified nucleotides, not a non-nucleotide spacer.'), ('UNG', 'Uracil-processing and PCR share one actual source bulk addition. Abasic fragments and other companions are not assumed physically removed by enzyme action alone.')))]
    content.append(panel(inventory('block').rows(),caption='Whole source blocking oligo; five terminal U residues are 2prime-O-methyl RNA'))
    content.extend(['<h2>Read setup</h2><p>'+'MiSeq Reagent Kit v2, 300-cycle kit, or NextSeq 75SE as reported. The MiSeq kit capacity does not assign an actual read/index mask; selected sample-index values and whole primer landing remain unassembled.'+'</p>','<h2>Source limits</h2><p>'+'The biological input boundary is isolated exosomes plus prepared antibody–DNA conjugates, not a complete exosome-isolation/conjugate-manufacturing graph. The selected exosome route does not use the starred STV-only adapter variants in SI1 Table1. Main primer names are not silently mapped to unverified table members. Unknown occupied circles, reporter extensions and read geometry remain explicit.'+'</p>','<details><summary>Full workflow, stock fractions and unresolved chemistry</summary>',WorkflowSection('Selected published preparation'),'</details></main>'])
    return assemble_page(g,content)

if __name__=='__main__':selftest()
