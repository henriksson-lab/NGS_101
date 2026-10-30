"""Published NJU-seq selected HeLa TR/NC preparation; private incomplete graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
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
    clean('fragment_precipitate','Source10ul3MNaOAc5.2/24ul5ugul glycogen/372ul96percentcoldethanol;literalamountsnotrepaired')
    allocated={a:g.product(a+'_rna','Source1ugRNAallocationfor'+a+'tail/MgRcomparison',K.MOLECULAR)for a in ('tr','nc')};rem=g.product('fragment_stock_remainder','Untransferredfragmentstock;actualamountunassigned',K.MOLECULAR);g.react('allocate_pair','SourceTR/NCpreparationsfromfragmentedRNAstock',S.SPLIT,((x,R.TEMPLATE),),(*allocated.values(),rem));g.mark_terminal(Terminal(rem.id,Disposition.STORED,'Untransferredstocknotassertedconsumed'))
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
        pnk=g.input(key('pnk'),'SourceT4PNK3phosphataseminus1ul/ATP10mM5ul/10xbuffer5ul/RNA39ul',K.CATALYST);y=g.product(key('phosphorylated'),'Source5phosphorylatedRNAclasses;3endphosphatesnotconvertedby3phosphataseminusassumption',K.MOLECULAR);g.react(key('phosphorylate'),'Source5phosphorylationAFTERMgR/controlcleanup',S.PHOSPHORYLATION,((start,R.TEMPLATE),(pnk,R.CATALYST)),(y,));start=y
        kit=g.input(key('kit'),'SourceNEBNextSmallRNALibraryPrepSetforIllumina;exacthistoricalversion/adapter/RT/PCRpartsunassembled',K.MOLECULAR);y=g.product(key('kit_product'),'SourcekitpreparedNGSmixture;fullintermediatesandwholelibraryunassembled',K.MOLECULAR);g.react(key('smallrna_kit'),'SourcehistoricalNEBsmallRNAkitboundary',S.MAGIC,((start,R.TEMPLATE),(kit,R.STRUCTURAL)),(y,),reason='Sourcekitpreparationpositive;wholemodifiedoligos,exactendstages/quantitativeallocationsandsealedproductsunassembled.');start=y
        cleanup('gel','Source6percentPAGE/fragmentslessthan50ntaswritten;insert-versus-wholelibrarysizeassignmentunassembled')
        lib=g.product(key('library'),'Source'+arm+'NGSlibrary;wholeends/primers/readmaskunassembled',K.LIBRARY);g.react(key('identity'),'Sourcefinalsequencinglibraryidentity',S.GENERIC,((start,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceIlluminaX-Ten;actualrun/read/indexcyclesunassigned. MgRtreatedTRandwithoutMgRNCareonepairedmethod,notextraprotocols.'))
    g.add_source_item(SourceItem('selected',SOURCE,'KnownpolyA-RNAtailorderingand5PNKsourcequalified;sealedend/readstatesunassembled.'));return g

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
    c('Privatepartialgraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
