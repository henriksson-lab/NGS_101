"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.'):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),reason=reason);return out

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

def protocol_graph():
    g=ProtocolGraph('RNA Tagging — associated-manuscript transcriptome-wide route',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source yeast expressing selected RBP–PUP2 fusion; target RNA and U-tag lengths unassigned',K.MOLECULAR,primary=True)
    x=op(g,'u_tag',x,'In-vivo RBP-associated PUP2 terminal uridylation; original transcripts retained',S.TAILING,catalyst=True,reagent='source RBP–PUP2 and cellular nucleotide pool',product='U-tagged and untagged RNA with all cell material; no perfect specificity or quantitative U-tag prediction')
    x=op(g,'lysis',x,'Source SDS/PCA/glass-bead lysis, ten vortex/ice cycles',S.EXTRACTION,reagent='RNA ISO buffer, acid-washed beads and PCA')
    x=clean(g,'organic',x,'Source bead removal, aqueous PCA then chloroform extraction and ethanol recovery')
    x=op(g,'dnase',x,'Source TURBO DNase1h37°C',S.DIGESTION,catalyst=True,reagent='source TURBO DNase',product='RNA plus DNA digestion products/reagents; DNA is not removed merely by digestion')
    x=clean(g,'genejet_rna',x,'Source GeneJetRNA purification')
    x=op(g,'polya',x,'Source Dynabeads mRNA polyA selection, approximately75µg totalRNA input',reagent='historical Dynabeads mRNA purification kit, edition unreviewed')
    x=op(g,'rrna',x,'Source yeast RiboZeroGold depletion after polyA selection',reagent='historical RiboZeroGold yeast kit, internal mechanism/remnant partition unreviewed')
    x=op(g,'gi_tail',x,'Source yeastPAP1200U/GTP0.5mM/ITP0.15mM90min37°C then another1200U30min37°C',S.TAILING,catalyst=True,reagent='yeast polyA polymerase plus explicit GTP and ITP; products contain riboinosine')
    x=clean(g,'tail_cleanup',x,'Source PCA/chloroform aqueous extraction and sodium-acetate/GlycoBlue ethanol recovery')
    primer=g.input('u_select','Associated-AM U-select RT primer Methods137; source trimming string separately preserved in note',K.OLIGO,state=oligo_scene('u_select',U_SELECT))
    enzyme=g.input('ssiii','Source SuperScriptIII50°C60min then85°C5min',K.CATALYST)
    cdna=g.product('first_dna','Source complementary firstDNA plus original RNA/enzymes; U-select primer contributes three A positions, not three independently observed U-tags',K.MOLECULAR)
    g.react('rt','Source U-select reverse transcription',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST)),(cdna,))
    cdna=op(g,'rnase_h',cdna,'Source RNaseH20min37°C',S.DIGESTION,catalyst=True,reagent='source RNaseH',product='cDNA plus RNA digestion fragments, undigested/carryover material; no automatic perfect RNA removal')
    cdna=clean(g,'cdna_cleanup',cdna,'Source GeneJetPCR purification, two32µl elutions')
    random=g.input('second_primer','Source adapter-tailed random hexamer; random priming is not a UMI',K.OLIGO,state=oligo_scene('second_primer',SECOND))
    klenow=g.input('klenow','Source exonuclease-minus Klenow37°C30min',K.CATALYST)
    second=g.product('second_dna','Source random-primed secondDNA products, original firstDNA and reaction remnants; priming boundaries unassigned',K.MOLECULAR)
    g.react('second_strand','Source Klenow second-strand synthesis',S.EXTENSION,((cdna,R.TEMPLATE),(random,R.STRUCTURAL),(klenow,R.CATALYST)),(second,))
    second=clean(g,'second_cleanup',second,'Source two1:1RNA CleanXP purifications')
    f=g.input('forward','Associated-AM RP1 primer',K.OLIGO,state=oligo_scene('rp1',RP1))
    r=g.input('reverse','Source sample-indexed PCR primer; actual sample index not assigned; no fabricated study value',K.OLIGO)
    taq=g.input('taq','Source GoTaqGreen cycling program Methods167',K.CATALYST)
    amp=g.product('amp','Source pooled PCR reactions including original templates/remnants; exact duplex/end/read geometry unassembled',K.MOLECULAR)
    g.react('pcr','Source GoTaqPCR, approximately9 parallel reactions pooled',S.PCR,((second,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(taq,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Random priming boundaries, selected index and exact full template assembly not assigned.'))
    amp=clean(g,'size',amp,'Source twice0.8:1RNA CleanXP size selection; synthetic-standard gel route not substituted')
    g.sequence(amp,SequencingPlan(acquisition='unspecified',limitation='AssociatedAMMethods181 explicitly HiSeq2500 PE50, R1gene identity/R2terminal tag; exact primer landing/end/index geometry and final-VOR applicability remain unassembled.'))
    g.add_source_item(SourceItem('pending','Main final-VOR and historical kit/predecessor/source-version reconciliation','Published main and named Methods119–181 vendor procedures; no claimed absence.'))
    return g

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
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'rna-tagging__10.1038+nmeth.3651'/'nmeth.3651_PMC4707952.xml.txt',"python3 tools/get_sources.py --slugs rna-tagging__10.1038+nmeth.3651; inspected original selected chemistry sections")
    if have(c,source,label="local original transcription guard"):
        c('source transcribes operational U-select oligo',U_SELECT in source.path.read_text(),True)
        c('source transcribes random second-strand primer',SECOND in source.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
