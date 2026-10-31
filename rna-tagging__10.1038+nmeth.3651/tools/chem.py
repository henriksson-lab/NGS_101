"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,molecular=False,additional=(),product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.'):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.MOLECULAR if molecular else K.CONDITION)
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    ports=((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.STRUCTURAL if molecular else R.CONDITION))+tuple(additional)
    if kind is S.EXTRACTION and molecular:
        condition=g.input(key+'_conditions',label+'; physical handling, not another reagent dose',K.CONDITION)
        ports+=((condition,R.CONDITION),)
    g.react(key,label,kind,ports,(out,),reason=reason);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc

from illumina import TRUSEQ_SMALL_RNA_RA5, TRUSEQ_SMALL_RNA_RA3, P5, P7
# Associated AM Methods137,163,167: original primer value, no trimming-string repair.
U_SELECT='G'+revcomp(TRUSEQ_SMALL_RNA_RA3)+'C'*9+'AAA'
SECOND=TRUSEQ_SMALL_RNA_RA5+'N'*6
RP1=P5+'GTTCAGAGTTCTACAGTCCGA'
SOURCE='Lapointe et al.2015 DOI10.1038/nmeth.3651; associated author manuscript PMC4707952 Methods119–181; final-VOR procedural applicability unverified.'
@dataclass(frozen=True)
class TaggedRNA:
    original:str
    u_count:int
    gi_tail:str=''
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGU') or type(self.u_count)is not int or self.u_count<0 or set(self.gi_tail)-set('GI'):raise ValueError('explicit native RNA, nonnegative U-tag and RNA G/I residues required')
    @property
    def residues(self):return self.original+'U'*self.u_count+self.gi_tail
    def tail(self,gi):
        if not gi or set(gi)-set('GI'):raise ValueError('supplied G/I tail must contain actual G or inosine residues')
        return replace(self,gi_tail=self.gi_tail+gi)
    def u_select_first_dna(self):
        # Explicit selected register only: source primer C9 pairs supplied G/I9,
        # source primer AAA covers U3. Those twelve positions are primer-derived,
        # never reverse-transcribed residues or independent observations of U-tags.
        if self.u_count<3 or len(self.gi_tail)!=9:
            raise ValueError('this selected U-select register requires at least3 terminalU and exactly9 supplied G/I residues; no automatic priming register prediction')
        upstream=self.original+'U'*(self.u_count-3)
        return U_SELECT+''.join({'A':'T','C':'G','G':'C','U':'A'}[b] for b in upstream[::-1])

def reverse_primer(index):
    if len(index)!=6 or set(index)-set('ACGT'):raise ValueError('source sample index is an explicitly supplied six-base DNA value')
    sc=Scene();sc.strand('sample_index_primer',[Segment('p7',P7),Segment('index',index,feature=feature('rna_tagging_sample_index','sample_index','fixed')),Segment('body','GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA')]);return sc

def graph_for_reaction_count(reaction_count=9):
    if type(reaction_count)is not int or not 1<=reaction_count<=10:
        raise ValueError('Illustrative PCR count must be positive and fit 5uL allocations from the source50uL eluate')
    g=ProtocolGraph('RNA Tagging — associated-manuscript transcriptome-wide route',authority=SOURCE,representation='illustrative_example',example_description='Representative PCR allocation count, default9, illustrates the source report of approximately9 reactions. This is not an asserted exact experimental count; sample index remains unassigned.')
    x=g.input('cells','Source yeast expressing selected RBP–PUP2 fusion; target RNA and U-tag lengths unassigned',K.BIOLOGICAL,primary=True)
    x=op(g,'u_tag',x,'In-vivo RBP-associated PUP2 terminal uridylation; original transcripts retained',S.MAGIC,reagent='Source in-vivo culture/expression conditions; PUP2 and nucleotide pool are already inside input cells, not fresh reagent doses',product='U-tagged and untagged RNA with all cell material; no perfect specificity or quantitative U-tag prediction')
    x=op(g,'lysis',x,'Source SDS/PCA/glass-bead lysis, ten vortex/ice cycles',S.EXTRACTION,reagent='RNA ISO buffer, acid-washed beads and PCA',molecular=True)
    x=clean(g,'organic',x,'Source bead removal, aqueous PCA then chloroform extraction and ethanol recovery')
    x=op(g,'dnase',x,'Source TURBO DNase1h37°C',S.DIGESTION,catalyst=True,reagent='source TURBO DNase',product='RNA plus DNA digestion products/reagents; DNA is not removed merely by digestion')
    x=clean(g,'genejet_rna',x,'Source GeneJetRNA purification')
    x=op(g,'polya',x,'Source Dynabeads mRNA polyA selection, approximately75µg totalRNA input',reagent='historical Dynabeads mRNA purification kit, edition unreviewed')
    x=op(g,'rrna',x,'Source yeast RiboZeroGold depletion after polyA selection',reagent='historical RiboZeroGold yeast kit, internal mechanism/remnant partition unreviewed')
    gt=g.input('gtp','Actual source0.5mM GTP initial tailing contribution',K.MOLECULAR)
    it=g.input('itp','Actual source0.15mM ITP initial tailing contribution; riboinosine retained',K.MOLECULAR)
    x=op(g,'gi_tail',x,'Initial G/I tailing90min37°C',S.TAILING,catalyst=True,reagent='Actual initial1200U yeastPAP',additional=((gt,R.STRUCTURAL),(it,R.STRUCTURAL)))
    x=op(g,'gi_tail_second',x,'Actual additional yeastPAP dose, another30min37°C; retained nucleotide mixture',S.TAILING,catalyst=True,reagent='Actual fresh second1200U yeastPAP dose; no fresh GTP/ITP dose asserted')
    x=clean(g,'tail_cleanup',x,'Source PCA/chloroform aqueous extraction and sodium-acetate/GlycoBlue ethanol recovery')
    primer=g.input('u_select','Associated-AM U-select RT primer Methods137; source trimming string separately preserved in note',K.OLIGO,state=oligo_scene('u_select',U_SELECT))
    enzyme=g.input('ssiii','Source SuperScriptIII50°C60min then85°C5min',K.CATALYST)
    rt_nts=g.input('rt_dntps','Actual source1uL10mM dNTP input for selective RT',K.MOLECULAR)
    cdna=g.product('first_dna','Source complementary firstDNA plus original RNA/enzymes; U-select primer contributes three A positions, not three independently observed U-tags',K.MOLECULAR)
    g.react('rt','Source U-select reverse transcription',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST),(rt_nts,R.STRUCTURAL)),(cdna,))
    cdna=op(g,'rnase_h',cdna,'Source RNaseH20min37°C',S.DIGESTION,catalyst=True,reagent='source RNaseH',product='cDNA plus RNA digestion fragments, undigested/carryover material; no automatic perfect RNA removal')
    cdna=clean(g,'cdna_cleanup',cdna,'Source GeneJetPCR purification, two32µl elutions')
    random=g.input('second_primer','Source adapter-tailed random hexamer; random priming is not a UMI',K.OLIGO,state=oligo_scene('second_primer',SECOND))
    klenow=g.input('klenow','Source exonuclease-minus Klenow37°C30min',K.CATALYST)
    nts=g.input('second_dntps','Actual source5uL10mM dNTP contribution to second-strand synthesis',K.MOLECULAR)
    second=g.product('second_dna','Source random-primed secondDNA products, original firstDNA and reaction remnants; priming boundaries unassigned',K.MOLECULAR)
    g.react('second_strand','Source Klenow second-strand synthesis',S.EXTENSION,((cdna,R.TEMPLATE),(random,R.STRUCTURAL),(klenow,R.CATALYST),(nts,R.STRUCTURAL)),(second,))
    second=clean(g,'second_cleanup',second,'Source two1:1RNA CleanXP purifications')
    aliquots=[g.product('pcr_'+str(i)+'_template','Illustrative reaction '+str(i)+': actual5uL source-template allocation',K.MOLECULAR)for i in range(1,reaction_count+1)]
    rest=[]
    if reaction_count*5<50:
        leftover=g.product('template_remainder','Unallocated template stock: '+str(50-5*reaction_count)+'uL in this illustrative allocation; later fate unassigned',K.MOLECULAR)
        rest.append(leftover)
    g.react('pcr_allocation','ILLUSTRATIVE — independent5uL template allocations from50uL eluate; source reports approximately9 reactions',S.SPLIT,((second,R.TEMPLATE),),tuple(aliquots+rest))
    for leftover in rest:
        g.mark_terminal(Terminal(leftover.id,Disposition.UNRESOLVED,'No source experiment-specific residual volume or later use asserted.'))
    products=[]
    for i,template in enumerate(aliquots,1):
        key='pcr_'+str(i)
        f=g.input(key+'_forward','Actual0.67uL10uM RP1 aliquot for this reaction',K.OLIGO,state=oligo_scene('rp1',RP1))
        r=g.input(key+'_reverse','Actual0.67uL10uM indexed-primer aliquot; same sample index unassigned, not a new sample',K.OLIGO)
        mix=g.input(key+'_mix','Actual8.33uL2xGoTaqGreen master-mix aliquot',K.MOLECULAR)
        water=g.input(key+'_water','Actual source2uL water for this PCR reaction',K.MOLECULAR)
        taq=g.product(key+'_taq','Polymerase component of this actual master-mix aliquot, not an extra dose',K.CATALYST)
        nucleotides=g.product(key+'_dntps','Nucleotide component of same master-mix aliquot; exact composition unassembled',K.MOLECULAR)
        other=g.product(key+'_buffer','Other copresent master-mix components; composition unassembled',K.MOLECULAR)
        g.react(key+'_components','Expose copresent components of one source-named master-mix aliquot',S.GENERIC,((mix,R.STRUCTURAL),),(taq,nucleotides,other),changes_molecule=False)
        product=g.product(key+'_product','Independent PCR product, original template and reagent remnants; full geometry unassembled',K.MOLECULAR)
        g.react(key,'Source GoTaq PCR of independent allocated template',S.PCR,((template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(taq,R.CATALYST),(nucleotides,R.STRUCTURAL),(other,R.STRUCTURAL),(water,R.STRUCTURAL)),(product,),geometry=PCRGeometry.unavailable(Scene(),reason='Random priming boundaries, sample index and full template assembly unassigned.'))
        products.append(product)
    amp=g.product('amp','Physically pooled independently amplified reactions and remnants; exact library geometry unassembled',K.MOLECULAR)
    g.react('pcr_pool','Source pools individual PCR reactions before cleanup',S.POOL,tuple((p,R.TEMPLATE)for p in products),(amp,))
    amp=clean(g,'size',amp,'Source twice0.8:1RNA CleanXP size selection; synthetic-standard gel route not substituted')
    g.sequence(amp,SequencingPlan(acquisition='unspecified',limitation='AssociatedAMMethods181 explicitly HiSeq2500 PE50, R1gene identity/R2terminal tag; exact primer landing/end/index geometry and final-VOR applicability remain unassembled.'))
    g.declare_reviewed_scope(('cells','u_tag_reagent','u_tag','gtp','itp','gi_tail_reagent','gi_tail','gi_tail_second_reagent','gi_tail_second','u_select','ssiii','rt_dntps','rt','rnase_h_reagent','rnase_h','second_primer','klenow','second_dntps','second_strand','pcr_allocation','pcr_pool','size_condition','size'),SOURCE)
    for i in range(1,reaction_count+1):
        key='pcr_'+str(i)
        g.declare_reviewed_scope(tuple(key+suffix for suffix in ('_forward','_reverse','_mix','_water','_components','')),SOURCE+'; actual per-reaction recipe with explicitly illustrative array count.')
    g.add_source_item(SourceItem('pending','Main final-VOR and historical kit/predecessor/source-version reconciliation','Published main and named Methods119–181 vendor procedures; no claimed absence.'))
    return g

def protocol_graph():
    return graph_for_reaction_count()

def selftest():
    from checks import Check
    c=Check();rna=TaggedRNA('ACGU',5).tail('GIGGGIGGG')
    c('U-tag and riboinosine remain physical RNA additions',rna.residues,'ACGU'+'U'*5+'GIGGGIGGG')
    c('selected U-select product retains primer-covered tail and copies upstreamRNA',rna.u_select_first_dna(),U_SELECT+'AAACGT')
    c('tailing retains the original transcript',rna.original,'ACGU')
    try:TaggedRNA('ACGU',1).tail('T');bad=False
    except ValueError:bad=True
    c('G/I tail cannot silently accept DNA T residues',bad,True)
    c('sample-index identity structured',bool(reverse_primer('ACGTAC').strands['sample_index_primer'].segs[1].feature),True)
    c('incomplete source/library remains unpublishable',bool(protocol_graph().publication_problems()),True)
    g=protocol_graph()
    c('Pool receives independent PCR products from allocated templates',all(g.producer_of(p.material).kind is S.PCR and g.producer_of(g.producer_of(p.material).inputs[0].material).id=='pcr_allocation' for p in g.steps['pcr_pool'].inputs),True)
    c('Whole-stock allocation creates no zero-volume remainder','template_remainder' in graph_for_reaction_count(10).materials,False)
    c.raises('Cannot allocate more than source50uL into5uL PCR templates',lambda:graph_for_reaction_count(11))
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'rna-tagging__10.1038+nmeth.3651'/'nmeth.3651_PMC4707952.xml.txt',"python3 tools/get_sources.py --slugs rna-tagging__10.1038+nmeth.3651; inspected original selected chemistry sections")
    if have(c,source,label="local original transcription guard"):
        c('source transcribes operational U-select oligo',U_SELECT in source.path.read_text(),True)
        c('source transcribes random second-strand primer',SECOND in source.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
