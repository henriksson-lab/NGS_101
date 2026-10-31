"""m6ACE selected protection event; source-limited graph, not a complete library."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Segment, Scene, feature, revcomp, panel
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_views import assemble_page, WorkflowSection
from page import head,info
PAGE_FILE="m6ace.html"
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, Terminal, Disposition, SequencingPlan
SOURCE='Koh et al.2019 DOI10.1038/s41467-019-13561-z; VOR Results39, Methods170 and Supplementary Data1'

@dataclass(frozen=True)
class TerminalC:
    base: str='C'
    sugar: str='2,3-dideoxyribose'
    three_prime_oh: bool=False
    def __post_init__(self):
        if (self.base,self.sugar,self.three_prime_oh)!=('C','2,3-dideoxyribose',False):
            raise ValueError('source terminal ddC is a blocked cytidine nucleotide')

@dataclass(frozen=True)
class Adapter:
    body: str=il.TRUSEQ_SMALL_RNA_RA3
    terminal: TerminalC=TerminalC()
    activation5: str='rApp'
    @property
    def template(self): return self.body+self.terminal.base
    def rt_primer(self): return revcomp(self.template)
    def construct(self):
        return Construct([Segment('adapter body',self.body),Segment('terminal ddC nucleotide',self.terminal.base)])

@dataclass(frozen=True)
class ProtectionEvent:
    """Conditional single-molecule path, not enzyme efficiency or raw spike preparation."""
    rna: str
    modified_position: int
    phosphate5: bool
    crosslinked: bool
    adapter: Adapter=Adapter()
    digested: bool=False
    washed: bool=False
    def __post_init__(self):
        if set(self.rna)-set('ACGU') or not 0<=self.modified_position<len(self.rna) or self.rna[self.modified_position]!='A':
            raise ValueError('explicit selected m6A must be an in-range RNA adenosine')
        if self.digested and not(self.phosphate5 and self.crosslinked):
            raise ValueError('selected XRN1 protection path requires 5′p and crosslinked antibody')
        if self.washed and not self.digested:raise ValueError('wash cannot precede this digest')
    def digest(self):
        if not(self.phosphate5 and self.crosslinked):raise ValueError('XRN1-ready protected substrate required')
        return replace(self,digested=True)
    @property
    def retained_rna(self):return self.rna[self.modified_position:] if self.digested else self.rna
    @property
    def soluble_nucleotides(self):return self.rna[:self.modified_position] if self.digested and not self.washed else ''
    def wash(self):
        if not self.digested:raise ValueError('digest first')
        return replace(self,washed=True)

def five_prime_adapter():
    return Construct([Segment('RA5 RNA body',il.TRUSEQ_SMALL_RNA_RA5.replace('T','U')),Segment('constrained eight-base UMI','NNNNNNNB',placeholder=True,feature=feature('m6ace_umi','umi','degenerate'))])

def primer_scene():
    a=Adapter();sc=Scene();sc.strand('adapter',[Segment('adapter including terminal ddC',a.template)],mod5='rApp',mod3='ddC: no 3′OH')
    sc.anneal('primer',[Segment('RT primer',a.rt_primer())],to='adapter',pair=('RT primer','adapter including terminal ddC'))
    return sc

SOURCE_OLIGO_EXCEPTIONS={'adapter3':'Source /3ddC/ is an appended physical blocked cytidine nucleotide. Generic source parsing omits this nucleotide; the model retains it. Whole modified original-order fixture and computed full RT-primer pairing verify the narrow transcription-format mapping. No source absence or chemical exception.'}
def unknown_scene():
    sc=Scene();sc.strand('context',[Segment('occupied product unassembled','[unassembled molecule]',placeholder=True,inferred=True)]);return sc

class Preparation:
    def __init__(self,g,x,prefix=''):self.g=g;self.x=x;self.prefix=prefix
    def key(self,k):return self.prefix+k
    def material(self,k,label,kind=K.MOLECULAR):
        m=self.g.input(self.key(k),label,kind);return(m,R.SELECTION if kind is K.SELECTION else R.STRUCTURAL)
    def cat(self,k,label):return(self.g.input(self.key(k),label,K.CATALYST),R.CATALYST)
    def cond(self,k,label):return(self.g.input(self.key(k),label,K.CONDITION),R.CONDITION)
    def op(self,k,label,kind,extra=(),chemical_state=None,geometry=None):
        out=self.g.product(self.key(k)+'_product','Actual '+label+' material; full occupied molecules/recovery/carryover unassembled',K.MOLECULAR,chemical_state=chemical_state)
        kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind is S.GENERIC:kw['changes_molecule']=False
        self.g.react(self.key(k),label,kind,((self.x,R.TEMPLATE),)+tuple(extra),(out,),**kw);self.x=out
    def clean(self,k,label,extra=()):
        out=self.g.product(self.key(k)+'_kept','Actual retained '+label+' fraction; recovery unassigned',K.MOLECULAR);waste=self.g.product(self.key(k)+'_waste','Actual excluded '+label+' workup fractions',K.WASTE)
        self.g.react(self.key(k),label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(extra),(out,waste));self.g.discard(waste,'Physical source workup reject; not a guarantee of purity');self.x=out

def pcr_supplies(b):
    g=b.g;k=b.key('pcr_supplies');stock=g.input(k,'Actual PCR aggregate supplies; source does not assign stock packaging/formulation here',K.MOLECULAR)
    enzyme=g.product(k+'_enzyme','Actual DNA polymerase functional component',K.CATALYST);nt=g.product(k+'_nts','Copresent nucleotide component of actual PCR supplies',K.MOLECULAR);buf=g.product(k+'_buffer','Copresent PCR buffer component; exact formulation unassembled',K.CONDITION)
    g.react(k+'_inventory','Copresent actual PCR components; not physical separation or a claimed master-mix factory',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nt,buf),changes_molecule=False)
    return((enzyme,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION))

def adapter3_scene():
    a=Adapter();sc=Scene();sc.strand('adapter',list(a.construct()),mod5='rApp',mod3='ddC: no3OH',label='Free source DNA adapter; terminal ddC is an additional physical cytidine')
    sc.ligation_end('adapter','adapter body',kind='single-strand',polarity='5-prime',owner='adapter',text='Known source5rApp activated donor; activationAMP released on intended ligation, ddC not an activation nucleotide')
    return sc

def protocol_graph():
    g=ProtocolGraph('m6ACE-seq — actual shared UV preparation and input/enriched libraries',authority=SOURCE,representation='published_instance')
    rna=g.input('rna','Actual source totalRNA; actual sequences/cohort identities unassembled',K.MOLECULAR,primary=True);b=Preparation(g,rna)
    b.clean('initial_ethanol','Actual initial ethanol precipitation to remove residual salt',(b.material('initial_ethanol_reagent','Actual ethanol precipitation materials'),))
    b.clean('polya','Actual Poly(A)Purist Mag selection',(b.material('polya_matrix','Actual source polyA selection materials',K.SELECTION),))
    b.clean('polya_ethanol','Actual ethanol recovery of polyA RNA',(b.material('polya_ethanol_reagent','Actual precipitation materials'),))
    b.op('fragment','RNA fragmentation70°C7.5min;120–150ntdistribution',S.FRAGMENTATION,(b.cond('fragment_buffer','Actual source RNA fragmentation buffer; composition unassembled'),))
    b.clean('fragment_ethanol','Actual fragment ethanol precipitation',(b.material('fragment_ethanol_reagent','Actual precipitation materials'),))
    b.op('pnk','Actual 10 U PNK37°C30min beforeATP',S.END_REPAIR,(b.cat('pnk_enzyme','Actual 10 U PNK'),))
    b.op('pnk_atp','Then actual1mMATP andadditional37°C30min',S.END_REPAIR,((g.materials['pnk_enzyme'],R.CATALYST),b.material('atp','Actual source1mMATP addition; carriedPNK notfreshenzyme'),),chemical_state={'bulkRNA':'Source end-treatedbulkRNA; intendedRNA5P/3OHclasses, occupancyunassembled','spike':'Not yet added; no earlierPNKassignedto laterpreligatedspike'})
    b.clean('end_cleanup','Oligo Clean and Concentrator recovery',(b.material('end_column','Actual source cleanup materials',K.SELECTION),))
    a=g.input('adapter3','Actual 200 pmol sourceDNA3adapter with5rApp/terminalddC',K.OLIGO,state=adapter3_scene(),chemical_state={'raw_order':'/5rApp/'+Adapter().body+'/3ddC/','terminal':'Physical2,3dideoxycytidinenucleotide; no3OH','limits':'Knownfreeadapter, notwhole occupied RNA or maturePCRcopy'})
    b.op('ligate3','Source3ligation25°C2h inATP-freebuffer',S.LIGATION,((a,R.STRUCTURAL),b.cat('rnl2','Actual 400 U truncatedT4RNA ligase2'),b.material('lig3_bsa','Actual60µg/mlBSAcarrier'),b.material('lig3_peg','Actual 12.5% PEG8000'),b.cond('lig3_buffer','Actual50mMTris7.5/10mMMg/10mMDTT ATP-freebuffer')),geometry=AdapterLigationGeometry.known('single-strand',adapter3_scene(),statement='Knownactivated5rAppdonor; source-intendedend-treatedRNA3OHacceptor. WholeoccupiedRNA/sealedjunction/enddistributionunassembled.'))
    b.clean('lig3_cleanup','Actual AMPure recovery of3ligatedbulkRNA',(b.material('lig3_beads','Actual AMPure cleanup materials',K.SELECTION),))
    spike=g.input('spike','Actual200pgalready3ligatedmethylatedRNAspike addedAFTERbulkPNK/ligationcleanup; free5endunassigned',K.MOLECULAR,chemical_state={'source_RNA':'UCUCCUUAGUAGCUCCUAAGm6AUCGUCGAGUUACACGACGACAUUGUUCCGACUGACG','adapter':'Sourcealready3ligatedstock; fullend/preparationunassembled','5end':'Unassigned; no priorbulkPNKretroactivelyapplied','XRN1':'Onlyeligible5P/crosslinkedclassfollows selectedprotectionevent, notallspikemolecules'})
    b.op('spike_add','Add200pgpreligatedspike to1–5µgligatedpolyARNAs',S.POOL,((spike,R.TEMPLATE),))
    b.op('denature','65°C5min thenice2min',S.THERMAL_TREATMENT,(b.cond('denature_condition','Actual RNA denaturation/ice schedule'),))
    b.op('antibody','Incubate with8µganti-m6A overnight4°C',S.BINDING,(b.material('antibody_protein','Actual selectedsourceanti-m6Aantibody; specificity/crosslinkefficiencynotassumed'),b.material('ip_inhibitor','Actual 1 U /µlRNasinPlus'),b.cond('ip_buffer','Actual150mMLiCl/10mMTris7.4/0.1%IGEPAL IPbuffer')))
    one=g.product('uv_one','Actual one 50 µl pre-UV aliquot',K.MOLECULAR);other=g.product('uv_other','Actual other pre-UV aliquot collection; individual50µlallocations/count/lastportionunassembled',K.MOLECULAR)
    g.react('uv_allocate','PhysicalpreUV50µlaliquotarray:oneclassandotherdistinctaliquotclasses',S.SPLIT,((b.x,R.TEMPLATE),),(one,other))
    uv=g.input('uv_conditions','Actual254nmUV six0.15J/cm²exposuresonaliquots',K.CONDITION)
    one_uv=g.product('uv_one_product','Actual crosslinked first aliquot class; unreacted/backgroundfractionsunassigned',K.MOLECULAR);other_uv=g.product('uv_other_product','Actual crosslinked other aliquot collection; notclonesofoneproduct',K.MOLECULAR)
    g.react('uv_crosslink','ActualUVtreatmentonphysicallydistinctaliquotarray; nofixedaliquotcount',S.CROSS_LINKING,((one,R.TEMPLATE),(other,R.TEMPLATE),(uv,R.CONDITION)),(one_uv,other_uv))
    recombined=g.product('uv_recombined','Actual recombined UV-treated RNA/antibody/spike mixture',K.MOLECULAR)
    g.react('uv_pool','Recombineactualtreatedaliquots AFTER UV',S.POOL,((one_uv,R.TEMPLATE),(other_uv,R.TEMPLATE)),(recombined,))
    input_rna=g.product('input_rna','Actual 1% post-UV input RNA control',K.MOLECULAR);enriched=g.product('enriched_rna','Actual 99% post-UV m6ACE fraction before bead capture',K.MOLECULAR)
    g.react('input_split','Actual post-UV recombined stock:1%input/99%m6ACE',S.SPLIT,((recombined,R.TEMPLATE),),(input_rna,enriched))
    beads=g.input('proteina_beads','Actual 1.2 mg DynabeadsProteinA blockedovernight4°CinIPbuffer',K.MOLECULAR);bsa=g.input('block_bsa','Actual 0.5 mg/ml BSA blockingcontribution',K.MOLECULAR)
    blocked=g.product('blocked_beads','Actual decanted BSA-blocked Protein A bead stock',K.SELECTION);block_waste=g.product('block_waste','Actual decanted blocking fluid',K.WASTE)
    g.react('block_beads','Actual parallel bead blocking and decant before capture',S.CLEANUP,((beads,R.TEMPLATE),(bsa,R.STRUCTURAL)),(blocked,block_waste));g.discard(block_waste,'Actual decanted blocking fluid')
    for prefix,start in (('input_',input_rna),('enriched_',enriched)):
        c=Preparation(g,start,prefix)
        if prefix=='enriched_':
            c.clean('capture','ProteinAcapture1.5h4°C andactualsequentialcoldwashes',((blocked,R.SELECTION),c.material('capture_washes','Actual250µl each washbuffer1/2/3/TE/Tris; molecularsalts/detergents/buffercomponents')))
            c.op('denature','DenaturebeadRNAin10µlTris65°C5min/ice2min',S.THERMAL_TREATMENT,(c.cond('denature_condition','Actual source preXRN1denaturation'),))
            c.op('xrn1','Source1UXRN1 37°C1h/shaking; selectedprotectedclassonly',S.DIGESTION,(c.cat('xrn1_enzyme','Actual 1 U XRN1'),c.material('xrn1_inhibitor','Actual 1 U /µlRNasinPlus'),c.cond('xrn1_buffer','Actual100mMLiCl/45mMTris8/10mMMg/1mMDTT')),chemical_state={'selectedbulkRNA':'Conditional5P/crosslinkedprotectedRNA startsatselectedm6A; crosslinkremnants/efficiencyunassigned','spike':'Starting5endunassigned; noforcedeligibility','solubleproducts':'Upstreamnucleotidesremaincopresentuntilphysicalwash'})
            c.clean('xrn1_wash','ActualpostXRN1washbuffer1/2/3/TE/Tris',(c.material('xrn1_washes','Actual source washsalts/detergents/buffers'),))
        c.op('proteinase','Actual PK/SDS elution50°C1.5hshaking',S.DEPROTEINIZATION,(c.cat('pk','Actual 1 mg/ml ProteinaseK'),c.material('elution_sds','Actual 1% SDS molecularcontribution'),c.cond('elution_buffer','Actual200mMNaCl/25mMTris8/2mMEDTA')),chemical_state={'RNA':'SourceinputorenrichedRNAwithpossiblecrosslink/peptideremnants; perfectdeproteinationnotasserted'})
        c.clean('ethanol','Actual RNA ethanol precipitation',(c.material('ethanol_reagent','Actual source ethanolprecipitationmaterials'),))
        five=g.input(prefix+'adapter5','Actual 5 pmol sourceRNA5adapterN7B',K.OLIGO,state=five_prime_adapter(),chemical_state={'polymer':'SourceRNAadapter','meaning':'ConstrainedeightntUMI NNNNNNNB, notN8','freeends':'Noextraadapterterminalmodificationsassigned'})
        c.op('ligate5','INFERRED — partialendview; sourceRNA5ligation16°C16h',S.LIGATION,((five,R.STRUCTURAL),c.cat('rnl1','Actual 10 U T4RNA ligase'),c.material('lig5_peg','Actual 12.5% PEG8000'),c.material('lig5_inhibitor','Actual 2 U/µl RNasin Plus')),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'SourceactualRNA5adapterligation; whole occupied RNA/freeadapterend/junctiongeometryunassembled. Source5Pbulkclassandspike unknown 5′ end retained.'))
        c.clean('lig5_cleanup','ActualOligoCleanConcentratorrecovery',(c.material('lig5_column','Actual source RNAcleanupmaterials',K.SELECTION),))
        primer=g.input(prefix+'rt_primer','Actual 5 pmol sourceRTprimer; anneal72°C2min/ice2min',K.OLIGO,state=Construct([Segment('computedsourceRTprimer',Adapter().rt_primer())]))
        rt_supply=g.input(prefix+'rt_supplies','Actual RT aggregate supplies; nucleotide/buffer amounts and stock packaging unassigned',K.MOLECULAR)
        rt_nt=g.product(prefix+'rt_nts','Copresent RT nucleotide functional class; no extra dose',K.MOLECULAR)
        rt_buf=g.product(prefix+'rt_buffer','Copresent RT buffer class; exact source formulation unassembled',K.CONDITION)
        g.react(prefix+'rt_inventory','Actual RT supply components; bookkeeping, not physical separation or a claimed premix',S.GENERIC,((rt_supply,R.TEMPLATE),),(rt_nt,rt_buf),changes_molecule=False)
        c.op('rt','Actual SSIII50°C1h; RNAtemplateretainedunlesssourceworkupchangesit',S.REVERSE_TRANSCRIPTION,((primer,R.STRUCTURAL),c.cat('ssiii','Actual 200 U SuperScriptIII'),(rt_nt,R.STRUCTURAL),(rt_buf,R.CONDITION)),chemical_state={'DNA':'ActualcDNApopulation; fullcopies/endlandingunassembled','RNA':'RTsourceRNAtemplatestillcopresent; noseparatevalidationRNaseborrowed','ddC':'OriginaladaptertemplateblockedddCretained; no3blockerassignedtocopiedDNA'})
        c.op('rt_stop','Actual70°C15minstop',S.THERMAL_TREATMENT,(c.cond('rt_stop_condition','Actual source 70°C15min'),))
        stock=g.input(prefix+'phusion','ActualPhusionHighFidelityPCRmastermixstockallocation',K.MOLECULAR);enzyme=g.product(prefix+'phusion_enzyme','CopresentPhusionpolymerasecomponent',K.CATALYST);nts=g.product(prefix+'phusion_nts','CopresentdNTPcomponentofsameactualmastermix',K.MOLECULAR);buffer=g.product(prefix+'phusion_buffer','CopresentPCRbuffercomponentofsameactualmastermix',K.CONDITION)
        g.react(prefix+'phusion_inventory','SamephysicalPhusionmastermixcomponents; noextraindependentdoses',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nts,buffer),changes_molecule=False)
        forward=g.input(prefix+'pcr_forward','Actual historical TruSeq forward primer functional class; whole order unassembled',K.OLIGO);reverse=g.input(prefix+'pcr_reverse','Actual historical TruSeq reverse primer functional class; whole order/index member unassembled',K.OLIGO)
        c.op('pcr','Actual 14-cycle library PCR',S.PCR,((enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='SourceactualPhusion/TruSeqPCR; wholehistoricalprimerorders/occupiedtemplate/copy/landingunassembled, no modernkitdefault.'))
        c.clean('ampure','Actual AMPure primer/adaptor-dimer removal',(c.material('ampure_beads','Actual source AMPurepurificationbeads',K.SELECTION),))
        load=g.product(prefix+'loading','Actual NextSeq PE75 loading contribution; amountunassigned',K.LIBRARY);rest=g.product(prefix+'library_stock','Untransferred prepared library stock; amount/fate unassigned',K.LIBRARY)
        g.react(prefix+'loading_allocate','Allocateactualsequencingcontribution',S.SPLIT,((c.x,R.TEMPLATE),),(load,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source prepared library stock; noinventedQCaliquot'))
        g.sequence(load,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':75,'Read 2':75},limitation='SourceNextSeqPE75; Read1first8ntN7BUMIprocessingexplicit. Fullreadprimerlanding/indexcycles/wholemoleculesunassembled.'))
    g.declare_reviewed_scope(COMMON_SCOPE,SOURCE)
    for prefix in ('input_','enriched_'):g.declare_reviewed_scope(tuple(prefix+k for k in BRANCH_SCOPE),SOURCE)
    g.declare_reviewed_scope(tuple('enriched_'+k for k in ENRICHED_SCOPE),SOURCE)
    return g.require_publishable()

COMMON_SCOPE=tuple("""rna initial_ethanol initial_ethanol_kept initial_ethanol_waste initial_ethanol_reagent polya polya_kept polya_waste polya_matrix polya_ethanol polya_ethanol_kept polya_ethanol_waste polya_ethanol_reagent fragment fragment_product fragment_buffer fragment_ethanol fragment_ethanol_kept fragment_ethanol_waste fragment_ethanol_reagent pnk pnk_product pnk_enzyme pnk_atp pnk_atp_product atp end_cleanup end_cleanup_kept end_cleanup_waste end_column adapter3 ligate3 ligate3_product rnl2 lig3_bsa lig3_peg lig3_buffer lig3_cleanup lig3_cleanup_kept lig3_cleanup_waste lig3_beads spike spike_add spike_add_product denature denature_product denature_condition antibody antibody_product antibody_protein ip_inhibitor ip_buffer uv_allocate uv_one uv_other uv_conditions uv_crosslink uv_one_product uv_other_product uv_pool uv_recombined input_split input_rna enriched_rna proteina_beads block_bsa blocked_beads block_waste block_beads""".split())
BRANCH_SCOPE=tuple("""proteinase proteinase_product pk elution_sds elution_buffer ethanol ethanol_kept ethanol_waste ethanol_reagent adapter5 ligate5 ligate5_product rnl1 lig5_peg lig5_inhibitor lig5_cleanup lig5_cleanup_kept lig5_cleanup_waste lig5_column rt_primer ssiii rt_supplies rt_inventory rt_nts rt_buffer rt rt_product rt_stop rt_stop_product rt_stop_condition phusion phusion_inventory phusion_enzyme phusion_nts phusion_buffer pcr_forward pcr_reverse pcr pcr_product ampure ampure_kept ampure_waste ampure_beads loading_allocate loading library_stock""".split())
ENRICHED_SCOPE=tuple("""capture capture_kept capture_waste capture_washes denature denature_product denature_condition xrn1 xrn1_product xrn1_enzyme xrn1_inhibitor xrn1_buffer xrn1_wash xrn1_wash_kept xrn1_wash_waste xrn1_washes""".split())

def render_page():
    g=protocol_graph();five=Scene();five.strand('RNA5',list(five_prime_adapter()),label='FreeRNA5adapter; sourceN7B constrains8ntUMI')
    return assemble_page(g,(head('m6ACE-seq — protected RNA to read-start libraries'),'<main class="wrap" data-protocol-complete="true"><h1>m6ACE-seq</h1><p>Antibody/UV protection and XRN1 enrich RNA fragments whose starts report protected modification sites. The retained RNA nucleotide remains present; this is a read-start signal, not a forced DNA mutation or deletion.</p><p><a href="https://doi.org/10.1038/s41467-019-13561-z">Koh et al. (2019)</a> · <a href="01_m6ace.html">Research note</a></p>',info('Published source-qualified preparation. INFERRED — selected 5′-phosphate/crosslinked protection events and intended ligation end classes do not establish complete digestion, crosslink efficiency, whole transcript sequence or occupied library geometry. The separately added preligated spike has an unassigned starting 5′ end; earlier bulk PNK is not silently applied to it.'),'<h2>Known free adapter components</h2>',panel(adapter3_scene().rows(),caption='Source activated DNA3 adapter with physical terminal ddC nucleotide; the blocker is not a C3 linker and is not inherited as a blocker by PCR copies.'),panel(primer_scene().rows(),caption='Computed full source RT-primer duplex includes the terminal ddC nucleotide as a template position paired to G; this free-component view is not a whole RNA library.'),panel(five.rows(),caption='Source RNA5 adapter with constrained NNNNNNNB UMI, not unrestricted N8.'),'<h2>Control allocation and reads</h2><p>Actual 50 µl UV aliquots are recombined before a 1% input/99% enriched split. Only the enriched branch receives Protein A capture/XRN1. Both branches continue through proteinase recovery, RNA5 ligation, SSIII, fourteen-cycle Phusion PCR, AMPure and NextSeq PE75. Read1 begins with the source eight-base UMI. Unknown full primer bodies, acquired index cycles and final landing remain explicit; no final QC aliquot or physical pool between controls is invented.</p>','<details><summary>Full source preparation, actual UV array and two DNA-library paths</summary>',WorkflowSection('Source preparation and physical fates'),'</details></main>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=Adapter();primer_scene()
    c('terminal ddC remains a template nucleotide',a.template.endswith('GC'),True)
    c('RT primer pairs terminal ddC and preserves blocked template',(a.rt_primer()[0],a.terminal.three_prime_oh),('G',False))
    p=ProtectionEvent('UACGAU',1,True,True).digest()
    c('digest retains selected modified A and downstream RNA',(p.retained_rna,p.soluble_nucleotides),('ACGAU','U'))
    c('nucleotide mass retained before wash',p.soluble_nucleotides+p.retained_rna,p.rna)
    c('wash physically separates digestion products',(p.wash().retained_rna,p.wash().soluble_nucleotides),('ACGAU',''))
    try:ProtectionEvent('UAC',1,False,True).digest();bad=False
    except ValueError:bad=True
    c('5′OH cannot silently enter selected XRN1 path',bad,True)
    c('Known donor retains physical blockedddC',protocol_graph().materials['adapter3'].chemical_state['terminal'],'Physical2,3dideoxycytidinenucleotide; no3OH')
    original=Source(sources_root()/'m6ace-seq__10.1038+s41467-019-13561-z'/'s41467-019-13561-z_supp_41467_2019_13561_MOESM4_ESM.xlsx.txt','Supplementary Data1')
    if have(c,original,label='source oligo workbook'):
        text=''.join(original.path.read_text().split())
        c('source adapter includes actual terminal ddC','/5rApp/'+a.body+'/3ddC/' in text,True)
        c('source RT primer matches computed full template',a.rt_primer() in text,True)
        c('source UMI preserves constrained B',(il.TRUSEQ_SMALL_RNA_RA5+'NNNNNNNB').replace('T','U') in text,True)
    g=protocol_graph()
    c('UVarray exists before crosslink',tuple(i.material for i in g.steps['uv_crosslink'].inputs if i.role is R.TEMPLATE),('uv_one','uv_other'))
    c('Pool combines distinct produced UVclasses',tuple(i.material for i in g.steps['uv_pool'].inputs if i.role is R.TEMPLATE),('uv_one_product','uv_other_product'))
    c('Actual postUV input/enriched split',tuple(g.steps['input_split'].outputs),('input_rna','enriched_rna'))
    c('Input doesnot inherit XRN1','input_xrn1' in g.steps,False)
    c('Spike doesnot inherit prior bulk PNK','Unassigned' in g.materials['spike'].chemical_state['5end'],True)
    c('Two actual DNA-library outputs',len(g.sequencing),2)
    c('Normal publication gate',g.publication_problems(),[])
    c.report(exit_on_fail=True)
if __name__=='__main__':selftest()
