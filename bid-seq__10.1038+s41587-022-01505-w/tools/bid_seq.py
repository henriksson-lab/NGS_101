"""BID-seq selected RNA-adduct/cDNA-deletion path, not a deletion-rate predictor."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Segment,Scene,feature,revcomp,panel
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_views import assemble_page,WorkflowSection
from page import head,info
PAGE_FILE="bid-seq.html"
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan
SOURCE='Dai et al.DOI10.1038/s41587-022-01505-w; published PMC10017504 Methods233–237, Results87–103 and analysis281'
RT_PRIMER=il.TRIM_SEEN_IN_READ1[:21]
RT_PRIMER=revcomp(RT_PRIMER)

@dataclass(frozen=True)
class SitePath:
    rna: str
    site: int
    adduct: bool=False
    skipped: bool|None=None
    rna_retained: bool=True
    def __post_init__(self):
        if set(self.rna)-set('ACGU') or not 0<=self.site<len(self.rna) or self.rna[self.site]!='U':raise ValueError('RNA-letter U at explicitly designated Ψ position required')
        if self.skipped and not self.adduct:raise ValueError('selected Ψ-BS deletion path requires adduct; backgrounds are outside this path')
    def treat(self):
        if self.skipped is not None:raise ValueError('treatment precedes copying')
        return replace(self,adduct=True)
    def copy(self,*,skip):
        if self.skipped is not None:raise ValueError('copy once')
        return replace(self,skipped=skip)
    @property
    def first_strand(self):
        if self.skipped is None:raise ValueError('copying has not occurred')
        template=self.rna[:self.site]+self.rna[self.site+1:] if self.skipped else self.rna
        return revcomp(template.replace('U','T'))

def adapter3():
    return Construct([Segment('3′ random window','N'*5,placeholder=True,feature=feature('bid_random_3','inline_barcode','random',note='Random source window; UMI interpretation remains inferred, not supplied by N5 alone')),
      Segment('inline barcode','ATCACG',feature=feature('bid_inline','inline_barcode','fixed')),
      Segment('3′ adapter DNA',il.TRIM_SEEN_IN_READ1[:21])])
def adapter5():
    return Construct([Segment('5′ RNA body',il.TRUSEQ_SMALL_RNA_RA5.replace('T','U')),Segment('5′ random window','N'*5,placeholder=True,feature=feature('bid_random_5','inline_barcode','random',note='Opposite-end random window, distinct from RNA3 adapter window; UMI semantics not asserted'))])
def adapter3_scene():
    sc=Scene();sc.strand('adapter',list(adapter3()),mod5='rApp',mod3='SpC3',label='Free source activated DNA adapter; C3 is a non-nucleotide blocker')
    sc.ligation_end('adapter','3′ random window',kind='single-strand',polarity='5-prime',owner='adapter',text='Source activated5rApp donor; activation AMP released during intended ligation, not part of DNA backbone')
    return sc

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

def protocol_graph():
    g=ProtocolGraph('BID-seq — actual shared preparation and input/treated libraries',authority=SOURCE,representation='published_instance')
    rna=g.input('rna','Actual selected biological-sample poly(A)+ RNA,10–20ng; independent biological replicate/cohort preparations not collapsed into this aliquot',K.MOLECULAR,primary=True)
    b=Preparation(g,rna)
    b.op('fragment','Fragment70°C14min',S.FRAGMENTATION,(b.cond('fragment_buffer','Actual source RNA fragmentation reagent'),))
    b.clean('fragment_cleanup','RNA Clean and Concentrator recovery',(b.material('fragment_column','Actual RNA cleanup materials',K.SELECTION),))
    b.op('end_repair','PNK first45min37°C in30µl',S.END_REPAIR,(b.cat('pnk_first','Actual 3 µlPNK'),b.cond('pnk_buffer','Actual 3 µl10xPNK buffer and water to30µl')))
    b.op('phosphorylate','Then freshPNK/ATP45min37°C',S.END_REPAIR,(b.cat('pnk_second','Actual fresh1.5µlPNK'),b.material('pnk_atp','Actual 1.5µl10 mM ATP')),chemical_state={'intended_source_ends':'Source3prime repair and5prime phosphorylation; actual occupancy/whole transcript/end distribution unassembled'})
    b.clean('end_cleanup','RNA cleanup; actual10µl eluate',(b.material('end_column','Actual RNA cleanup materials',K.SELECTION),))
    donor=g.input('rna3_adapter','Actual 1 µl20µM source RNA3 adapter; free DNA stock pattern',K.OLIGO,state=adapter3_scene(),chemical_state={'source_ends':'5rApp activation and non-nucleotide3SpC3 block','source_index':'PrintedATCACG inline barcode; full final library/index architecture unassembled'})
    b.op('rna3_ligation','RNA3 ligation after70°C2min/ice:25°C2h then16°C10h',S.LIGATION,((donor,R.STRUCTURAL),b.cat('rnl2','Actual 2 µlRnl2truncatedKQ'),b.material('rnl2_peg','Actual 7.5µl50% PEG8000'),b.material('rnl2_inhibitor','Actual 1 µlSUPERaseIn'),b.cond('rnl2_buffer','Actual 2.5µl10xligasebuffer;10µlRNA+1µladapter')),geometry=AdapterLigationGeometry.known('single-strand',adapter3_scene(),statement='Source activated5rApp donor; source-intended repairedRNA3OH. Whole occupied acceptor, releasedAMP distribution and sealed junction unassembled.'))
    b.op('deadenylate','Dilute reaction to47µl; add2µldeadenylase30°C30min',S.DIGESTION,(b.cat('deadenylase','Actual 2 µl5prime deadenylase; AMP/byproducts remain copresent until physical cleanup'),))
    b.op('recjf_digest','Then1µlRecJf37°C30min',S.DIGESTION,(b.cat('recjf','Actual 1 µlRecJf'),))
    b.clean('rna3_cleanup','RNA cleanup; actual9.3µl eluate',(b.material('rna3_column','Actual RNA cleanup materials',K.SELECTION),))
    five=g.input('rna5_adapter','Actual 1.2µl10µM source RNA5 adapter; no undocumented terminal modification assigned',K.OLIGO,state=adapter5(),chemical_state={'polymer':'RNA source ordering includesU; exact free terminal modifications unassigned','meaning':'N5window distinct from opposite-end random window; no automatic UMI assumption'})
    b.op('rna5_ligation','INFERRED — partial end view; source RNA5 ligation after70°C2min/ice;25°C8h',S.LIGATION,((five,R.STRUCTURAL),b.cat('rnl1','Actual 1 µlT4 RNA ligase 1'),b.material('rnl1_atp','Actual 1 µl25 mM ATP'),b.material('rnl1_peg','Actual 10 µl50% PEG8000'),b.cond('rnl1_buffer','Actual 2.5µl10xligasebuffer')),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),), 'Source actualRNA5adapter ligation to end-treatedRNA; whole occupied acceptor/donor and sealed junction not assembled. Intended adapter3OH/sourceRNA5P class does not establish all substrate termini.'))
    b.clean('rna5_cleanup','RNA cleanup; actual10µl eluate',(b.material('rna5_column','Actual RNA cleanup materials',K.SELECTION),))
    input_rna=g.product('input_rna','Actual 1.5µl of10µl preparedRNA for input library',K.MOLECULAR);treated=g.product('treated_rna','Actual 8.5µl of10µl preparedRNA for treated library',K.MOLECULAR)
    g.react('source_split','Actual shared10µlRNA allocation:1.5µlinput and8.5µltreated',S.SPLIT,((b.x,R.TEMPLATE),),(input_rna,treated))
    for prefix,start in (('input_',input_rna),('treated_',treated)):
        c=Preparation(g,start,prefix)
        if prefix=='treated_':
            c.op('bisulfite','Add45µlfreshBIDBSreagent;70°C3h',S.CHEMICAL_CONVERSION,(c.material('bs_reagent','Actual fresh 2.4MNa2SO3/0.36MNaHSO3 reagent; molecular reactants'),c.cond('bs_conditions','Source70°C3h')),chemical_state={'RNA_selected_site':'Ψ-BS adduct retains RNA residue/backbone; conversion/degradation fractions unassigned'})
            c.clean('bs_capture','Column binding/firstwash',(c.material('bs_binding','Actual75µlwater/270µlRNAbindingbuffer/400µlethanol; column and200µlwash contributions'),))
            c.op('desulfonate','Actual200µlRNA DesulphonationBuffer oncolumn1h roomtemperature',S.CHEMICAL_CONVERSION,(c.material('desulfonation_reagent','Actual source column desulfonation reagent; composition unassembled'),c.cond('desulfonation_condition','Source1h roomtemperature')),chemical_state={'RNA_selected_site':'SourceΨ-BS selected adduct path, RNA not deleted'})
            c.clean('bs_cleanup','Two700µlwashes/elute10.5µl',(c.material('bs_wash','Actual source RNAwashbuffer andwater'),))
        else:c.op('dilute','Dilute actual1.5µlinput to10.5µl',S.GENERIC,(c.material('dilution_water','ActualRNase-free water addition to10.5µl'),))
        rp=g.input(prefix+'rt_primer','Actual 1 µl2µM SR RT primer;65°C2min/ice annealing',K.OLIGO,state=Construct([Segment('source RT primer',RT_PRIMER)]))
        c.op('rt','SSIV50°C1h after source primerannealing',S.REVERSE_TRANSCRIPTION,((rp,R.STRUCTURAL),c.cat('ssiv','Actual 1 µlSuperScriptIV'),c.material('rt_nts','Actual 2 µl10 mM dNTP'),c.material('rt_dtt','Actual 1 µl100mMDTT'),c.material('rt_inhibitor','Actual 0.5µlRNaseOUT'),c.cond('rt_buffer','Actual 4 µl5xSSIVbuffer')),chemical_state={'RNA_site':'Ψ-BS'if prefix=='treated_'else'Ψ','RNA_retained':True,'cDNA':'Treated selected adduct can yield a skipped complementary base or nondeletion copy; background/unconverted/motif fractions unassigned. RNA itself is not physically deleted.'})
        c.op('rnase_h','Add1µlRNaseH37°C20min',S.DIGESTION,(c.cat('rnaseh_enzyme','Actual 1 µlRNaseH'),),chemical_state={'cDNA':'Selected deletion/nondeletionDNA classes remain distinct from RNA residue chemistry','RNA_fragments':'Copresent until physical DNA-column workup'})
        c.op('heat','70°C5min',S.THERMAL_TREATMENT,(c.cond('heat_conditions','Actual source70°C5min'),))
        c.clean('dna_cleanup','DNA cleanup;20µlstockstored−80°C',(c.material('dna_column','Actual DNA Clean and Concentrator materials',K.SELECTION),))
        selected=g.product(prefix+'pcr_template','Actual 4 µl of20µlcDNAstock for one representative sourcePCR allocation',K.MOLECULAR);rest=g.product(prefix+'cdna_stock','Nominal16µluntransferred cDNA stock; actual further PCR count unassigned',K.MOLECULAR)
        g.react(prefix+'pcr_allocate','Actual 4 µl perPCR allocation from20µlstock; no invented replicate count',S.SPLIT,((c.x,R.TEMPLATE),),(selected,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'UntransferredsourcecDNAstock; no disposalorcloning'));c.x=selected
        forward=g.input(prefix+'pcr_forward','Actual sourceSR Primerfor Illumina functionalclass; whole historical orderunassembled',K.OLIGO);reverse=g.input(prefix+'pcr_reverse','Actual sourceNEBNextindexedprimer functional class; whole body/member unassembled',K.OLIGO,chemical_state={'sample_index':'Sourceindexedprimerclass, exactvalue/width/coordinate unassigned'})
        c.op('pcr','Actual 15-cycle production PCR',S.PCR,pcr_supplies(c)+((forward,R.STRUCTURAL),(reverse,R.STRUCTURAL)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='ActualhistoricalSR/NEBNextindexedPCR, wholeoccupiedcDNA/primerorders/landingunassembled. No vendorrecipe or fixedsourceenzymeformulation inferred.'))
        c.clean('gel','Actual 3.5% low-melting agarose library purification',(c.material('gel_material','Actual sourcegelpurificationmaterials',K.SELECTION),))
        load=g.product(prefix+'loading','Actual selected NovaSeq6000library loading contribution; amountunassigned',K.LIBRARY);libstock=g.product(prefix+'library_stock','Untransferredprepared library stock; amount/fate unassigned',K.LIBRARY)
        g.react(prefix+'loading_allocate','Allocate actual source sequencing contribution',S.SPLIT,((c.x,R.TEMPLATE),),(load,libstock));g.mark_terminal(Terminal(libstock.id,Disposition.UNRESOLVED,'Sourcepreparedlibrary remainder'))
        g.sequence(load,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':100},limitation='SourceNovaSeq6000SE100; wholeprimerlanding, acquiredindexcycles andstudyindexmembersunassembled. SelectedRNAadductisnotitselfaDNAdeletion.'))
    g.declare_reviewed_scope(COMMON_SCOPE,SOURCE)
    for prefix in ('input_','treated_'):g.declare_reviewed_scope(tuple(prefix+k for k in BRANCH_SCOPE),SOURCE)
    g.declare_reviewed_scope(tuple('treated_'+k for k in TREATED_SCOPE),SOURCE);g.declare_reviewed_scope(('input_dilute','input_dilute_product','input_dilution_water'),SOURCE)
    return g.require_publishable()

COMMON_SCOPE=tuple("""rna fragment fragment_product fragment_buffer fragment_cleanup fragment_cleanup_kept fragment_cleanup_waste fragment_column end_repair end_repair_product pnk_first pnk_buffer phosphorylate phosphorylate_product pnk_second pnk_atp end_cleanup end_cleanup_kept end_cleanup_waste end_column rna3_adapter rna3_ligation rna3_ligation_product rnl2 rnl2_peg rnl2_inhibitor rnl2_buffer deadenylate deadenylate_product deadenylase recjf_digest recjf_digest_product recjf rna3_cleanup rna3_cleanup_kept rna3_cleanup_waste rna3_column rna5_adapter rna5_ligation rna5_ligation_product rnl1 rnl1_atp rnl1_peg rnl1_buffer rna5_cleanup rna5_cleanup_kept rna5_cleanup_waste rna5_column source_split input_rna treated_rna""".split())
BRANCH_SCOPE=tuple("""rt_primer ssiv rt_nts rt_dtt rt_inhibitor rt_buffer rt rt_product rnase_h rnase_h_product rnaseh_enzyme heat heat_product heat_conditions dna_cleanup dna_cleanup_kept dna_cleanup_waste dna_column pcr_allocate pcr_template cdna_stock pcr_forward pcr_reverse pcr_supplies pcr_supplies_inventory pcr_supplies_enzyme pcr_supplies_nts pcr_supplies_buffer pcr pcr_product gel gel_kept gel_waste gel_material loading_allocate loading library_stock""".split())
TREATED_SCOPE=tuple("""bisulfite bisulfite_product bs_reagent bs_conditions bs_capture bs_capture_kept bs_capture_waste bs_binding desulfonate desulfonate_product desulfonation_reagent desulfonation_condition bs_cleanup bs_cleanup_kept bs_cleanup_waste bs_wash""".split())

def render_page():
    g=protocol_graph();sc=Scene();sc.strand('RNA5',list(adapter5()),label='Free source RNA5 adapter; terminal modifications not assigned from its name')
    return assemble_page(g,(head('BID-seq — RNA adduct and DNA readout'),'<main class="wrap" data-protocol-complete="true"><h1>BID-seq</h1><p>Optimized bisulfite chemistry creates a pseudouridine–bisulfite RNA adduct. SSIV can omit a complementary base while copying that adduct; the RNA residue is not deleted. Shared adapter-ligated RNA is physically divided into input and treated library branches.</p><p><a href="https://doi.org/10.1038/s41587-022-01505-w">Dai et al. (2022)</a> · <a href="01_bid_seq.html">Research note</a></p>',info('Published defining-source preparation, not a successor protocol. INFERRED — source-intended acceptor/end classes and selected adduct/deletion paths do not determine whole transcript sequences, conversion rates or exact occupied library molecules. Unknown historical primer bodies, full landing and index cycles remain unassembled.'),'<h2>Known free adapter stocks</h2>',panel(adapter3_scene().rows(),caption='Source activated DNA3 adapter with non-nucleotide C3 block; N5 random window and printed ATCACG inline barcode do not establish a whole final library.'),panel(sc.rows(),caption='Source RNA5 adapter with a distinct N5 window; molecule-UMI interpretation remains unassigned.'),'<h2>Allocations and read setup</h2><p>The actual 10 µl prepared RNA stock splits into 1.5 µl input and 8.5 µl treated. Both routes reach 10.5 µl before RT. Each 20 µl cDNA stock contributes 4 µl to a representative 15-cycle PCR allocation, leaving 16 µl nominal stock; additional PCR count is unassigned. Gel-purified libraries use NovaSeq6000 SE100. No extra QC aliquot or physical pool between controls is invented.</p>','<details><summary>Full source preparation, both libraries and physical fates</summary>',WorkflowSection('Shared preparation and distinct input/treated outputs'),'</details></main>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=SitePath('AGUGA',2).treat().copy(skip=True)
    c('RNA adduct survives a selected cDNA deletion',(p.rna,p.adduct,p.first_strand,p.rna_retained),('AGUGA',True,'TCCT',True))
    c('same adduct permits a nondeletion path',SitePath('AGUGA',2).treat().copy(skip=False).first_strand,'TCACT')
    try:SitePath('AGUGA',2).copy(skip=True);bad=False
    except ValueError:bad=True
    c('specific adduct-deletion relation cannot silently describe untreated RNA',bad,True)
    c('opposite-end identifiers retain distinct identities',len({s.feature.id for s in [*adapter3(),*adapter5()] if s.feature}),3)
    original=Source(sources_root()/'bid-seq__10.1038+s41587-022-01505-w'/'s41587-022-01505-w_PMC10017504.xml.txt','published Methods233–237')
    if have(c,original,label='source library oligos'):
        text=''.join(original.path.read_text().split())
        c('RT primer source transcription',RT_PRIMER in text,True)
        c('3′adapter activation and non-nucleotide block preserved','5′App-NNNNNATCACG'+il.TRIM_SEEN_IN_READ1[:21]+'-3SpC3' in text,True)
        c('5′RNA adapter random barcode source',(il.TRUSEQ_SMALL_RNA_RA5+'NNNNN').replace('T','U') in text,True)
    g=protocol_graph()
    c('Actual common preparedRNA split',tuple(g.steps['source_split'].outputs),('input_rna','treated_rna'))
    c('Distinct freshPNK dose provenance',g.materials['pnk_first'].id!=g.materials['pnk_second'].id,True)
    c('Treated bisulfite is a molecular reactant',g.materials['treated_bs_reagent'].kind,K.MOLECULAR)
    c('Input branch avoids bisulfite','input_bisulfite' in g.steps,False)
    c('Actual production template allocations',tuple(g.steps[k+'pcr'].inputs[0].material for k in ('input_','treated_')),('input_pcr_template','treated_pcr_template'))
    c('Both nominal cDNAstock remainders remain unresolved',tuple(g.terminals[k+'cdna_stock'].disposition for k in ('input_','treated_')),(Disposition.UNRESOLVED,Disposition.UNRESOLVED))
    c('Two actual DNA-library outputs',len(g.sequencing),2)
    c('Normal publication gate',g.publication_problems(),[])
    c.report(exit_on_fail=True)
if __name__=='__main__':selftest()
