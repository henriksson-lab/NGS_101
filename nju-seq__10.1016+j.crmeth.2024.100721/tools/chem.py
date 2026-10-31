"""Published NJU-seq selected HeLa TR/NC preparation; private incomplete graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition,SourceCoverage
from adapter_ligation import AdapterLigationGeometry
TITLE='NJU-seq: selected HeLa TR/NC pair'
SOURCE='PublishedPMC10985248 Methods504/508,Fig1;oneoriginalcappedpass. Nm-VAQ,analyticalsubstratesandotherassaysnotmerged;historicalNEBkitinternalsunassembled,noabsenceclaim.'
TAIL_WRITTEN='pAGCUAAAAAAAAAAAAp'
TAIL='AGCU'+'A'*12

def tail_inventory():
    sc=Scene();sc.strand('polyA_ssRNA',[Segment('source_RNA_tail',TAIL,note='LiteralpAGCUAAAAAAAAAAAAp;RNApolymerandbothterminalpmarksretained')],mod5='P',mod3='P');return sc

def end_inventory():
    sc=tail_inventory();sc.strand('fragment',[Segment('source_fragment','[fragmented RNA end unassembled]',placeholder=True)],label='Fragment3endstateunassembled;noassumed3OH/endhealing');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');x=g.input('hela_rna','SourceHeLatotalRNA15ug;othercases10ugnotmerged',K.MOLECULAR,primary=True)
    probe=g.input('rrna_probe','Source1ulH/M/RrRNAprobein29ulRNA;wholeDNAorderingunassembled',K.OLIGO)
    rnase=g.input('rnaseh','SourceRNaseH1ul/10xbuffer4ul/water4ulinitialreaction;repeatadditionsand5/3assignmentunassembled',K.CATALYST)
    condition=g.input('deplete_condition','Source95C2min→22C.1Csec/5min;37C30min;repeat5/3timesaswritten',K.CONDITION)
    y=g.product('depleted_mix','SourceRNA/rRNAcleavageproducts/probe/reagentmixture;depletionnotperfect',K.MOLECULAR)
    g.react('rrna_deplete','Sourcehybridization/RNaseHcycles;exact5/3branchassignmentunassembled',S.MAGIC,((x,R.TEMPLATE),(probe,R.STRUCTURAL),(rnase,R.CATALYST),(condition,R.CONDITION)),(y,),reason='Repeatedsourcehybridization/RNaseHpositive;iterationassignment/freshstockamounts/wholeprobesunassembled.');x=y
    dnase=g.input('dnase','Source2UDNaseI37C30min',K.CATALYST);y=g.product('probe_digested','SourceRNAanddigestedDNAprobecompanions',K.MOLECULAR);g.react('digest_probe','SourceDNaseIafterrRNAdepletion',S.DIGESTION,((x,R.TEMPLATE),(dnase,R.CATALYST)),(y,));x=y
    def clean(key,label):
        nonlocal x
        b=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','SourcephysicallyretainedRNApreparation;yield/carryoverunassigned',K.MOLECULAR);w=g.product(key+'_waste','Actualphysicalworkupfraction',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalfractionnotcompletechemicalselection');x=y
    clean('deplete_clean','SourceVAHTSRNACleanBeads')
    frag=g.input('fragment_stock','SourceThermoRNAfragmentationreagents95C10min',K.CONDITION);y=g.product('fragmented','SourcefragmentedRNA;endchemistry/lengthdistributionunassembled',K.MOLECULAR);g.react('fragment','SourceRNAfragmentation',S.FRAGMENTATION,((x,R.TEMPLATE),(frag,R.CONDITION)),(y,));x=y
    precipitation=g.input('fragment_precipitants','Actual10ul3MNaOAc5.2/24ul5ugul glycogen/372ul96percent cold ethanol; literal source quantities',K.MOLECULAR)
    precip_mix=g.product('fragment_precipitation_mix','Fragmented RNA plus source solvent/salt/carrier before physical recovery',K.MOLECULAR)
    g.react('add_precipitants','Add actual source precipitation reagents',S.GENERIC,((x,R.TEMPLATE),(precipitation,R.STRUCTURAL)),(precip_mix,),changes_molecule=False)
    x=precip_mix
    clean('fragment_precipitate','Source10ul3MNaOAc5.2/24ul5ugul glycogen/372ul96percentcoldethanol;literalamountsnotrepaired')
    allocated={a:g.product(a+'_rna','Source1ugRNAallocationfor'+a+'tail/MgRcomparison',K.MOLECULAR)for a in ('tr','nc')};rem=g.product('fragment_stock_remainder','Untransferredfragmentstock;actualamountunassigned',K.MOLECULAR);g.react('allocate_pair','SourceTR/NCpreparationsfromfragmentedRNAstock',S.SPLIT,((x,R.TEMPLATE),),(*allocated.values(),rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferred stock; amount and later fate unassigned, no storage instruction inferred'))
    for arm,start in allocated.items():
        def key(k):return arm+'_'+k
        def cleanup(k,label):
            nonlocal start
            b=g.input(key(k+'_condition'),label,K.CONDITION);y=g.product(key(k+'_kept'),'Sourcephysicallyretainedmaterial;purity/recoveryunassigned',K.MOLECULAR);w=g.product(key(k+'_waste'),'Actualsourcewash/gelunretainedfraction',K.WASTE);g.react(key(k),label,S.CLEANUP,((start,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalworkupnotperfectmolecularseparation');start=y
        tail=g.input(key('tail'),'Source3ul25uMpolyA-ssRNAadaptor',K.ADAPTER,state=tail_inventory(),chemical_state={'polymer':'RNA','five_prime':'P','three_prime':'P','written':TAIL_WRITTEN,'scope':'Bothliteralpmarksretained;completefragmentendhandling/MgRentrygeometryunassembled'})
        c=g.input(key('tail_ligation_stock'),'Source16Covernightligation;exactcatalystidentity/cofactorstockunassembled',K.CATALYST);y=g.product(key('tailed'),'Sourcetail-ligatedandunligatedRNAclasses;wholejunction/endprocessingunassembled',K.MOLECULAR)
        geom=AdapterLigationGeometry('unavailable','unavailable',(end_inventory(),),'SourcepolyA-ssRNAtail5P/3Pmarksprinted;fullfragment3end,ligationreactionandphosphatehandlingunassembled;notasourceabsenceclaim.')
        g.react(key('tail_ligate'),'SourcepolyA-ssRNAadaptorligationbeforeMgR',S.LIGATION,((start,R.TEMPLATE),(tail,R.STRUCTURAL),(c,R.CATALYST)),(y,),geometry=geom);start=y
        if arm=='tr':
            c=g.input(key('mgr'),'SourceMgR30ug/60ul20mMTris8/100mMKCl/.01mMZnCl2',K.CATALYST);kind=S.DIGESTION;role=R.CATALYST
        else:c=g.input(key('mock'),'SourcesameprocesswithoutMgRprotein;buffer-onlycomparison',K.CONDITION);kind=S.THERMAL_TREATMENT;role=R.CONDITION
        y=g.product(key('digest_product'),'Source'+arm+'RNApopulation;MgRpausingatNm+1notassertedperfectorwholeRNAidentity',K.MOLECULAR);g.react(key('digest'),'Source37C30minMgRorNC',kind,((start,R.TEMPLATE),(c,role)),(y,));start=y;cleanup('rna_clean','SourceRNACleanandConcentrator5')
        pnk=g.input(key('pnk'),'Actual1ulT4PNK3prime-phosphatase-minus; source39ulRNAinput',K.CATALYST)
        atp=g.input(key('atp'),'Actual5ul10mM ATP contribution',K.MOLECULAR)
        pnkbuf=g.input(key('pnk_buffer'),'Actual5ul10xT4PNK ReactionBuffer',K.CONDITION)
        y=g.product(key('phosphorylated'),'Source5phosphorylatedRNAclasses;3endphosphatesnotconvertedby3phosphataseminusassumption',K.MOLECULAR)
        g.react(key('phosphorylate'),'Source5phosphorylationAFTERMgR/controlcleanup',S.PHOSPHORYLATION,((start,R.TEMPLATE),(pnk,R.CATALYST),(atp,R.STRUCTURAL),(pnkbuf,R.CONDITION)),(y,));start=y
        kit=g.input(key('kit'),'SourceNEBNextSmallRNALibraryPrepSetforIllumina;exacthistoricalversion/adapter/RT/PCRpartsunassembled',K.MOLECULAR);y=g.product(key('kit_product'),'SourcekitpreparedNGSmixture;fullintermediatesandwholelibraryunassembled',K.MOLECULAR);g.react(key('smallrna_kit'),'SourcehistoricalNEBsmallRNAkitboundary',S.MAGIC,((start,R.TEMPLATE),(kit,R.STRUCTURAL)),(y,),reason='Sourcekitpreparationpositive;wholemodifiedoligos,exactendstages/quantitativeallocationsandsealedproductsunassembled.');start=y
        cleanup('gel','Source6percentPAGE/fragmentslessthan50ntaswritten;insert-versus-wholelibrarysizeassignmentunassembled')
        lib=g.product(key('library'),'Source'+arm+'NGSlibrary;wholeends/primers/readmaskunassembled',K.LIBRARY);g.react(key('identity'),'Sourcefinalsequencinglibraryidentity',S.GENERIC,((start,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceIlluminaX-Ten;actualrun/read/indexcyclesunassigned. MgRtreatedTRandwithoutMgRNCareonepairedmethod,notextraprotocols.'))
    g.add_source_item(SourceItem('selected',SOURCE,'KnownpolyA-RNAtailorderingand5PNKsourcequalified;sealedend/readstatesunassembled.'));return _publication_coverage(g)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();sc=tail_inventory()
    c('WholeliteralpolyAtail','p'+TAIL+'p','pAGCUAAAAAAAAAAAAp')
    c('Bothsourcepmarksretained',(sc.strands['polyA_ssRNA'].mod5,sc.strands['polyA_ssRNA'].mod3),('P','P'))
    c('TailprecedesMgR',g.steps['tr_digest'].inputs[0].material,'tr_tailed')
    c('NCcontainsnoMgRenzyme',[(i.material,i.role)for i in g.steps['nc_digest'].inputs if i.role==R.CATALYST],[])
    c('PNKfollowsactualMgRworkup',g.steps['tr_phosphorylate'].inputs[0].material,'tr_rna_clean_kept')
    c('PNKnotdephosphorylation',g.steps['tr_phosphorylate'].kind,S.PHOSPHORYLATION)
    c('TRandNCbothsequenced',sum(t.disposition==Disposition.SEQUENCING for t in g.terminals.values()),2)
    src=Source(sources_root()/'nju-seq__10.1016+j.crmeth.2024.100721'/'j.crmeth.2024.100721_PMC10985248.xml.txt','https://pmc-oa-opendata.s3.amazonaws.com/PMC10985248.1/PMC10985248.1.xml;Methods508')
    if have(c,src):c('Wholemodifiedtailagainstselectedoriginal',TAIL_WRITTEN in src.path.read_text(),True)
    c('Both source PNK reactions consume actual ATP',all(arm+'_atp' in {i.material for i in g.steps[arm+'_phosphorylate'].inputs} for arm in ('tr','nc')),True)
    c('Untransferred fragment stock has unassigned fate',g.terminals['fragment_stock_remainder'].disposition,Disposition.UNRESOLVED)
    c('Source-qualified publication gate',g.publication_problems(),[]);c.report()


PAGE_FILE='nju-seq.html'
SOURCE_TARGETS=(('Publishedmain504/508 selectedHeLatotalRNA/rRNAdepletion/DNase/fragmentation/actualprecipitantamounts/matchedTRNC1ugallocations. 5/3repeatassignment remains sourcequalified, not silentlynormalized.', ('hela_rna', 'rrna_probe', 'rnaseh', 'deplete_condition', 'rrna_deplete', 'dnase', 'digest_probe', 'deplete_clean_condition', 'deplete_clean', 'fragment_stock', 'fragment', 'fragment_precipitants', 'add_precipitants', 'fragment_precipitate_condition', 'fragment_precipitate', 'allocate_pair')), ('Samepublishedmain508 actualTRtailbeforeMgR/RNAClean/5primePNKATP/historicalNEBsmallRNAkit/PAGE/XTen. Bothliteralterminalpmarks retained; wholeends/kitintermediates/readmaskunassembled.', ('tr_tail', 'tr_tail_ligation_stock', 'tr_tail_ligate', 'tr_mgr', 'tr_digest', 'tr_rna_clean_condition', 'tr_rna_clean', 'tr_pnk', 'tr_atp', 'tr_pnk_buffer', 'tr_phosphorylate', 'tr_kit', 'tr_smallrna_kit', 'tr_gel_condition', 'tr_gel', 'tr_identity')), ('Samepublishedmain508 says exactsameprocess forNCwithoutMgRprotein. Separateactualcontrolproducer; no pooling or enzyme introducedintoNC.', ('nc_tail', 'nc_tail_ligation_stock', 'nc_tail_ligate', 'nc_mock', 'nc_digest', 'nc_rna_clean_condition', 'nc_rna_clean', 'nc_pnk', 'nc_atp', 'nc_pnk_buffer', 'nc_phosphorylate', 'nc_kit', 'nc_smallrna_kit', 'nc_gel_condition', 'nc_gel', 'nc_identity', 'nc_library')))

def _publication_coverage(g):
    g.cover_source_item(SourceCoverage('selected','modelled','tr_library'))
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
    content=[head('NJU-seq'),'<main class="wrap" data-protocol-complete="true"><h1>'+'NJU-seq'+'</h1><p>'+'MgR digestion profiles RNA endpoint patterns associated with endogenous 2′-O-methylation. A poly(A)-RNA tail is ligated before digestion, and matched treated and no-MgR preparations continue through phosphorylation and small-RNA library preparation.'+'</p><p><a href="'+'https://doi.org/10.1016/j.crmeth.2024.100721'+'">Defining published article</a> · <a href="01_chem.html">Research note and source scope</a></p>',caveat('<strong>Source-qualified preparation.</strong> Actual source-known inputs and process order are retained. Unknown occupied molecules, yields and whole library/read geometry remain explicit; local reagent views do not prove a unique reaction product.'),'<h2>Preparation</h2><ol>'+'<li>Start with the source 15 µg HeLa total RNA. Perform rRNA-probe/RNase H depletion and DNase treatment, then recover and fragment RNA.</li><li>Use the literal source salt/glycogen/ethanol precipitation recipe. Allocate separate 1 µg contributions to treated and negative-control preparations; any stock remainder has an unresolved fate.</li><li>Ligate the source poly(A)-ssRNA tail before MgR. Treat TR with MgR; prepare NC by the same process without MgR protein. Purify each independently.</li><li>Apply the source 3′-phosphatase-minus T4 PNK with ATP for 5′ phosphorylation. Continue through the historical NEBNext small-RNA kit interface and separate PAGE recovery before sequencing.</li>'+'</ol>','<h2>Known molecular parts</h2>',table(('Part','Role and source limit'),(('Poly(A)-ssRNA tail', 'Literal pAGCUAAAAAAAAAAAAp ordering, with both terminal phosphate marks preserved. The source fragment recipient ends and whole sealed junction remain unassembled.'), ('MgR', 'RNA digestion enzyme; exact pausing/occupancy is not assumed perfect. The RNA modification itself is not converted into a different base by the model.'), ('T4 PNK', 'Actual post-digestion 5′ phosphorylation with ATP. The 3′-phosphatase-minus enzyme is not used to silently heal all source RNA 3′ ends.'), ('Small-RNA library interface', 'Source historical NEBNext preparation is positive, while its full oligos, internal RT/PCR recipe and products remain unassembled. No default kit constructor is imported.')))]
    content.append(panel(tail_inventory().rows(),caption='Whole source poly(A)-ssRNA tail; both terminal p marks retained, recipient RNA ends unassembled'))
    content.extend(['<h2>Read setup</h2><p>'+'Illumina X-Ten; actual run/read/index cycle counts remain unassigned. DNA read positions report surviving RNA endpoints; they are not a direct base conversion or a proof that every RNA endpoint corresponds to a methylated residue.'+'</p>','<h2>Source limits</h2><p>'+'This page follows published main504/508. The source says depletion is repeated “5/3” times without a resolved branch assignment here; that expression is retained. The source PAGE text says fragments smaller than 50 nt; insert versus whole-library size interpretation remains unassembled. Native RNA tail end marks and unresolved kit geometry are retained, not repaired from another methylation method. Nm-VAQ and analytical substrate assays are separate source branches.'+'</p>','<details><summary>Full workflow, stock fractions and unresolved chemistry</summary>',WorkflowSection('Selected published preparation'),'</details></main>'])
    return assemble_page(g,content)

if __name__=='__main__':selftest()
