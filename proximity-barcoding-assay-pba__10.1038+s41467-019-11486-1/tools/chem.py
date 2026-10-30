"""Published PBA exosome preparation; private incomplete graph, no MPX inheritance."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
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
    t=g.input('template','Sourceligationtemplate',K.OLIGO,state=inventory('template'));lig=g.input('circle_ligase','SourceT4ligase.01Uul/phi29buffer/ATP1mM/100nMeach37C30min',K.CATALYST);circle=g.product('circle_mix','Sourcecircularizedpadlockplusligationtemplate/unligatedcompanions; closuregeometryunassembled',K.MOLECULAR)
    g.react('circle','Sourcepadlockcircularization',S.CIRCULARIZATION,((pad,R.TEMPLATE),(t,R.STRUCTURAL),(lig,R.CATALYST)),(circle,),geometry=CircularizationGeometry.known('splinted single-strand ligation',statement='Publishedmain167ligates5primephosphorylatedpadlockonsuppliedligationtemplate; sourceTable1wholeorderingsknown. Fullpairedclosure/productstateunassembled;noevent/yieldclaim.'))
    exo=g.input('exonucleases','SourceExoI.2Uul/ExoIII2Uul37C30minthen85C20min',K.CATALYST);dig=g.product('circle_enriched','Sourceexonuclease-treatedcirclemixture; actualsurvivingprimers/linearfragmentsunassigned',K.MOLECULAR);g.react('exo','Sourceunligatedpadlockdepletion',S.DIGESTION,((circle,R.TEMPLATE),(exo,R.CATALYST)),(dig,))
    adj=g.input('primer_adjustment','SourceRCAprimersameorderingasligationtemplate adjusted100nM;additionalamountunassigned',K.OLIGO,state=inventory('template'));b=g.input('reanneal_condition','Source37C20min',K.CONDITION);templ=g.product('rca_template','Sourcecirclecomponentafterreannealing',K.MOLECULAR);primer=g.product('retained_rca_primer','Co-retainedsourceRCAprimercomponent;sametube,notphysicalfraction',K.OLIGO,state=inventory('template'));g.react('reanneal','SourceRCAprimeradjustment/reannealing',S.HYBRIDIZATION,((dig,R.TEMPLATE),(adj,R.STRUCTURAL),(b,R.CONDITION)),(templ,primer))
    phi=g.input('phi29','Sourcephi29.1Uul/dAGCU1mM37C10min',K.CATALYST);rca=g.product('rca_reaction','Sourceuracil-containingRCAproductspluscompanions; repeats/geometryunassembled',K.MOLECULAR);g.react('rca','SourceRCA',S.EXTENSION,((templ,R.TEMPLATE),(primer,R.STRUCTURAL),(phi,R.CATALYST)),(rca,))
    h=g.input('rca_heat','Source65C10min',K.CONDITION);rcastock=g.product('rca_stock','Sourceheat-terminatedRCAstockkeptminus20C;totalvolumeunassigned',K.MOLECULAR);g.react('rca_stop','SourceRCAheatstop',S.THERMAL_TREATMENT,((rca,R.TEMPLATE),(h,R.CONDITION)),(rcastock,));take=g.product('rca_selected','Source25ulRCAproductat1nM;stockdilution/transferamountunassigned',K.MOLECULAR);rem=g.product('rca_remainder','UntransferredRCAstock;amountunassigned',K.MOLECULAR);g.react('rca_aliquot','SourceRCAstockallocation',S.SPLIT,((rcastock,R.TEMPLATE),),(take,rem));g.mark_terminal(Terminal(rem.id,Disposition.STORED,'Sourcekeptstock;untransferredfraction'))
    probes=g.input('probe_pool','Sourcepreconjugated40probe/38proteinantibodypool20nMeach;SMCC/reduction/dialysisupstreamunassembled',K.MOLECULAR)
    op('stain','Sourceantibody-exosomebinding4Covernight',S.BINDING,'SourcePBAbufferinclsalmonDNA/BSA','SourceexosomeswithboundandfreeAOCpool',((probes,R.STRUCTURAL),))
    selected=g.product('one_ul','Source1ulstainedstock',K.MOLECULAR);left=g.product('stain_remainder','Untransferredstainedstock;amountunassigned',K.MOLECULAR);g.react('stain_aliquot','Actualsource1ulallocation',S.SPLIT,((x,R.TEMPLATE),),(selected,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferredstocklateruseunassigned'));x=selected
    op('dilute','Source1ulcomplexesdilutedin25ulPBAbuffer aswritten',S.GENERIC,'SourcePBAbuffer','Sourcedilutedcomplexstock;25ultotalvsdiluentnotresolved')
    aliquot=g.product('capture_20ul','Source20ulcaptureinput',K.MOLECULAR);left=g.product('dilution_remainder','Sourceuntransferreddilution;volumeunassigned',K.MOLECULAR);g.react('capture_aliquot','Actualsource20ulallocation',S.SPLIT,((x,R.TEMPLATE),),(aliquot,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferreddilution'));x=aliquot
    ctb=g.input('ctb_wells','SourceSTVcoatedtubewith50ngbiotinCTBovernightthen2washes;preparedcapturestock',K.SELECTION);b=g.input('capture_condition','Source15minRT',K.CONDITION);y=g.product('captured','SourceCTBboundexosomes/AOCs;occupancy/purityunassigned',K.MOLECULAR);w=g.product('uncaptured','Physicalunretainedcapturematerial',K.WASTE);g.react('capture','Sourceexosomecapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(ctb,R.SELECTION),(b,R.CONDITION)),(y,w));g.discard(w,'Sourcecapturedpreparationnotallvesicles');x=y;wash('wash1','Source2x100ulwashingbuffer')
    block=g.input('block','Source500nMblockingoligo',K.OLIGO,state=inventory('block'),chemical_state={'terminal':'mUmUmUmUmU','scope':'2primeO-methylRNA bases;notordinarydU'})
    op('bind_rca','Source25ulRCA1nM/block500nM37C15min',S.HYBRIDIZATION,'Sourcebindingtemperature/buffer','SourceRCA-boundreporterpreparation',((take,R.STRUCTURAL),(block,R.STRUCTURAL)));wash('wash2','Source2x100ulwashingbuffer')
    # AOC primers are components already present on captured exosomes, not free additions.
    ext=g.input('t4polymerase','Source1.25UT4DNApolymerase/25ulbuffer/dNTP100uMRT15min',K.CATALYST);y=g.product('extended','SourceAOCextensioncopiescomplexTag; someunextendedoligosmayundergo3primeexonucleolysis',K.MOLECULAR);g.react('extend','SourceproximityextensiononboundRCAtemplate;fullgeometryunassembled',S.MAGIC,((x,R.TEMPLATE),(ext,R.CATALYST)),(y,),reason='SourceboundAOC/RCAcopyingpositive; exactsubstrate/retainedAOCprimercomponentandfullcopiedproductunassembled.');x=y;wash('wash3','SourcewashbeforePCR')
    # Explicit components of each single source bulk addition, not extra additions.
    def bulk(key,label,enzyme=False):
        stock=g.input(key+'_bulk',label,K.MOLECULAR)
        f=g.product(key+'_forward','Forwardprimercomponentofsamebulkstock;exactorderingunassembled',K.OLIGO)
        r=g.product(key+'_reverse','Reverseprimercomponentofsamebulkstock;exactorderingunassembled',K.OLIGO)
        b=g.product(key+'_buffer','Buffer/dNTP/dye/DMSOcomponentsofsamebulkstock',K.CONDITION)
        outs=[f,r,b];enz=None
        if enzyme:
            enz=g.product(key+'_enzyme','SourcePhusion.02Uul+UNG.02Uulcomponentsofsamebulkaddition',K.CATALYST);outs.append(enz)
        g.react(key+'_components','Representsamebulkstockcomponents;notphysicalseparation',S.GENERIC,((stock,R.TEMPLATE),),outs,changes_molecule=False)
        return f,r,b,enz
    f,r,b,enz=bulk('first_mix','Source30ulPhusionHF/Phusion.02Uul/UNG.02Uul/dNTP/SYBR/DMSO/pba-fwd+rev500nMeach;exactoligomappingunassembled',True)
    y=g.product('pcr_product','SourceUNG-treated/amplifiedextensionproductreactionincludingretainedenzymes;abasicRCAcompanionsnotphysicallyremoved',K.MOLECULAR)
    g.react('pcr','SourcecombinedUNG37C15min→95C2min→30PCRcycles95C15s60C30s72C30s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.CONDITION)),(y,),reason='SourceoneUNG/Phusionbulkadditionandorderedprogramme;exactprimer/productassignmentunassembled.');x=y
    take=g.product('index_1ul','Source1ulwholeamplificationreactionforindexPCR',K.MOLECULAR);left=g.product('pcr_remainder','Sourceuntransferredampliconreactionstock',K.MOLECULAR);g.react('index_aliquot','Actualsource1ulallocation',S.SPLIT,((x,R.TEMPLATE),),(take,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'UntransferredPCRstock'));x=take
    def retained(key):
        nonlocal x
        template=g.product(key+'_template','DNAsubstratecomponentofselectedwholePCRreaction;notaphysicalfraction',K.MOLECULAR)
        enz=g.product(key+'_retained_enzyme','Co-transferredpriorPhusion/UNGcomponent;amountandactivityunassigned,noexternaladdition',K.CATALYST)
        g.react(key+'_carryover_components','Representco-retainedreactioncomponents;notnewstockaddition',S.GENERIC,((x,R.TEMPLATE),),(template,enz),changes_molecule=False);x=template;return enz
    enz=retained('index')
    f,r,b,_=bulk('index_mix','Source10ulindexPCRmixcontaining10uMeachfwd-index/rev-index;newpolymeraseidentity/amountnotassigned')
    y=g.product('index_product','Source2cycleindexedPCRreactionincludingco-transferredenzymes;fullproductunassembled',K.MOLECULAR)
    g.react('index','Source1ulspikedin10ulindexmix/2cycles95C15s60C60s72C60s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.CONDITION)),(y,),reason='SourceindexPCRpositive;co-transferredpriorcatalystamount/activityandcompleteprimer/productsunassembled,newenzymeadditionnotinferred.');x=y
    take=g.product('library_selected','SelectedindexedPCRreactionfor20folddilution;transferamountunassigned',K.MOLECULAR);left=g.product('index_remainder','UntransferredindexedPCRreactionstock;amountunassigned',K.MOLECULAR);g.react('library_aliquot','Actualsourceindexedstockallocation',S.SPLIT,((x,R.TEMPLATE),),(take,left));g.mark_terminal(Terminal(left.id,Disposition.STORED,'Untransferredindexedstock'));x=take
    enz=retained('library')
    f,r,b,_=bulk('library_mix','Sourcefreshmix1xPhusionHFbuffer/dNTP/SYBR/DMSO/library-fwd+rev500nM;buffernameisnotafreshenzymeaddition')
    y=g.product('library_pcr_product','Source20folddilution/15cyclelibraryPCRreaction;fullproducts/catalystactivityunassembled',K.MOLECULAR)
    g.react('library_pcr','Source20folddilutionintonewmixthen15cycles95C15s60C30s72C30s',S.MAGIC,((x,R.TEMPLATE),(enz,R.CATALYST),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.CONDITION)),(y,),reason='SourcefreshmixlistsPhusionHFbufferandprimers;newpolymeraseadditionnotinferred.Co-transferredpriorenzymeshaveunassignedamount/activity.');x=y
    # Multiple sample pooling is not expanded into fabricated other samples.
    op('pool_boundary','SourcePCRproductpoolingwithotherstudyproductsunassembled',S.MAGIC,'Sourcepoolingprogramme','Selectedsourceampliconcontribution;otherpoolmembersunassembled');wash('purify','SourceQIAquickPCRpurification')
    lib=g.product('library','SourcePBAsequencinglibrary;wholeprimerlayoutunassembled',K.LIBRARY);g.react('library_identity','Physicallibraryidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMiSeqv2 300cyclesorNextSeq75SE;perrunexactread/indexmasksunassigned. Bioinformaticminimumreadfiltersarecomputational,notDNAremoval.'))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Physicalexosome/RCAtemplatecopyingknown;exactcircle/extension/readgeometryunassembled.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('WholepublishedAb1',AB1,'AAGCGTGTGATTCCTAGTAATGAACAATGGNNNNNNNNTGCTATTATGATGTCTCAGGT')
    c('Wholepublishedpadlock',PADLOCK,'GCTATTATGATGTCTCAGGUNNNNNNNNNNNNNNNAGAUCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTAGTGCTGGATGATCGTCC')
    c('Wholepublishedtemplate',TEMPLATE,'CTGAGACATCATAATAGCGGACGATCATCCAGCACT')
    c('WholeoriginalblockingRNA',BLOCK,'CTAGTGCTGGATGATCGTCCUUUUU')
    c('RCAconsumesproducedprimer',[(i.material,i.role)for i in g.steps['rca'].inputs if i.role in {R.STRUCTURAL,R.SELF_PRIMER}],[('retained_rca_primer',R.STRUCTURAL)])
    c('RCAstockremainderexplicit',g.terminals['rca_remainder'].disposition,Disposition.STORED)
    c('Captureusesproduced20ul',g.steps['capture'].inputs[0].material,'capture_20ul')
    c('Sourcebindthenwashbeforeextension',g.steps['extend'].inputs[0].material,'wash2_kept')
    c('UNGPhusionaddedonesourcemix',[(i.material,i.role)for i in g.steps['pcr'].inputs if i.role==R.CATALYST],[('first_mix_enzyme',R.CATALYST)])
    c('Indexconsumesretainednotnewcatalyst',[(i.material,i.role)for i in g.steps['index'].inputs if i.role==R.CATALYST],[('index_retained_enzyme',R.CATALYST)])
    c('Libraryconsumesretainednotbufferasnewenzyme',[(i.material,i.role)for i in g.steps['library_pcr'].inputs if i.role==R.CATALYST],[('library_retained_enzyme',R.CATALYST)])
    c('Firstbulkprimerscontributemolecules',[(i.material,i.role)for i in g.steps['pcr'].inputs if i.role==R.STRUCTURAL],[('first_mix_forward',R.STRUCTURAL),('first_mix_reverse',R.STRUCTURAL)])
    c('Indexedstockremainderexplicit',g.terminals['index_remainder'].disposition,Disposition.STORED)
    c('Source1ulindexingstockremainder',g.terminals['pcr_remainder'].disposition,Disposition.STORED)
    c('Physicalselectionnotreadfilter',g.steps['purify'].kind,S.CLEANUP)
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
