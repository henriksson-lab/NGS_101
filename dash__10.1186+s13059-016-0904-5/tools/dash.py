"""DASH selected duplex break and source-owned downstream RNA-seq library processing."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
from pcr import PCRGeometry
SOURCE='Gu et al.2016 DOI10.1186/s13059-016-0904-5; published PMC4778327 Fig1, Methods191–211; published Additionalfile1 guide table'
GUIDE='ATTTTCAGTGTATTGCTTTG'

@dataclass(frozen=True)
class Fragment:
    dna: str
    outer_sites: frozenset[str]
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or not isinstance(self.outer_sites,frozenset) or self.outer_sites-{'left','right'}:raise ValueError('explicit DNA and immutable declared outer-adapter sites required')
    @property
    def paired_adapter_eligible(self):return self.outer_sites==frozenset({'left','right'})

@dataclass(frozen=True)
class SelectedCas9Break:
    fragment: Fragment
    recognition20: str
    pam_start: int
    def __post_init__(self):
        if len(self.recognition20)!=20 or set(self.recognition20)-set('ACGT'):raise ValueError('explicit20base source recognition domain required')
        if type(self.pam_start) is not int or not 20<=self.pam_start<=len(self.fragment.dna)-3:raise ValueError('in-range PAM following complete recognition domain required')
        if self.fragment.dna[self.pam_start-20:self.pam_start]!=self.recognition20 or self.fragment.dna[self.pam_start+1:self.pam_start+3]!='GG':raise ValueError('this selected exact-match NGG path does not represent mismatched/non-NGG cleavage')
    @property
    def cut_position(self):return self.pam_start-3
    def products(self):
        i=self.cut_position;f=self.fragment
        return Fragment(f.dna[:i],f.outer_sites&{'left'}),Fragment(f.dna[i:],f.outer_sites&{'right'})

def protocol_graph():
    g=ProtocolGraph('DASH — published RNA-seq depletion branch',authority=SOURCE,representation='published_instance')
    dna=g.input('tagmented_cdna','Source HeLa/CSF NuGEN-amplified, Nextera-partial-adapter cDNA; historical full polymers unassembled',K.MOLECULAR,primary=True)
    cas9=g.input('cas9','Source recombinant SpCas9 protein before guide assembly',K.MOLECULAR)
    guides=g.input('guide_pool','Source54 mitochondrial-rRNA sgRNAs; full RNA-start/scaffold composition unassembled',K.OLIGO)
    buffer=g.input('rnp_buffer','Source Cas9 buffer; RNP assembly37°C10min',K.CONDITION)
    rnp=g.product('rnp','Source assembled Cas9/sgRNA RNP population',K.CATALYST)
    g.react('assemble_rnp','Source RNP assembly; no guide-pool stoichiometry inferred from one selected target',S.BINDING,((cas9,R.TEMPLATE),(guides,R.STRUCTURAL),(buffer,R.CONDITION)),(rnp,))
    digest=g.product('digest','Source37°C2h digestion mixture: cut and uncut DNA plus RNP; both cut fragments retained',K.MOLECULAR,chemical_state={'DNA_fate':'no physical discard at cleavage','adapter_relation':'each selected cut product retains one original outer site','other_products':'unassigned'})
    g.react('cut','Cas9 treatment of adapter-bearing library DNA; no total depletion/off-target prediction',S.DIGESTION,((dna,R.TEMPLATE),(rnp,R.CATALYST)),(digest,))
    pk=g.input('proteinase','Source ProteinaseK1µl,37°C15min',K.CATALYST)
    proteolyzed=g.product('proteolyzed','DNA and proteolysis products remain in reaction; no automatic protein-fragment disposal',K.MOLECULAR)
    g.react('proteolyze','Source rRNA-library Cas9 workup branch',S.DEPROTEINIZATION,((digest,R.TEMPLATE),(pk,R.CATALYST)),(proteolyzed,))
    def clean(key,x,label,procedure):
        reagent=g.input(key+'_reagent',procedure,K.SELECTION)
        keep=g.product(key+'_kept',label,K.MOLECULAR)
        waste=g.product(key+'_waste','Physical '+key+' rejection fractions; exact molecular recovery unquantified',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(keep,waste));g.discard(waste,'Actual physical separation, without assuming complete species-specific removal');return keep
    extracted=clean('organic_extract',proteolyzed,'Recovered aqueous DNA mixture','Three phenol/chloroform/isoamyl extractions then one chloroform extraction')
    recovered=clean('ethanol_recovery',extracted,'Source precipitated/resuspended DNA','Source sodiumacetate/linearacrylamide/ethanol;4°C45min pellet and70%ethanol wash')
    f=g.input('index_forward','Source dual-index forward primer; exact historical polymer/index member unreviewed',K.OLIGO)
    r=g.input('index_reverse','Source dual-index reverse primer; exact historical polymer/index member unreviewed',K.OLIGO)
    kapa=g.input('kapa','Source KAPA library-amplification chemistry',K.CATALYST)
    def pcr(key,x,label):
        out=g.product(key+'_product',label,K.MOLECULAR)
        g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(kapa,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(Scene(),reason='Historical complete partial-adapter/primer/strand geometry and index members unassembled.'));return out
    amplified=pcr('index_pcr',recovered,'Source10–12cycle dual-index amplification; cut/original templates still present, no molecular deletion by PCR')
    selected=clean('size_selection',amplified,'Source360–540bp DNA fraction','BluePippin360–540bp selection; not an all-cut-fragments-rejected claim')
    purified=clean('dna_cleanup',selected,'Source Zymo-cleaned DNA','DNA Clean and Concentrator-5')
    reamplified=pcr('q_pcr',purified,'Source qPCR-guided reamplification, stopped in exponential phase; exact cycle count not assigned')
    other=g.input('other_libraries','Other source separately prepared, quantified libraries for sequencing pool',K.LIBRARY)
    pooled=g.product('pooled_library','Pooled source DASH/control libraries; mixing ratios and residual material unassigned',K.LIBRARY)
    g.react('pool','Source library pooling/quantification handoff',S.POOL,((reamplified,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,))
    g.sequence(pooled,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':135,'Read 2':135},limitation='Source HiSeq4000PE135; exact final-library/run-primer/index landing remains unassembled.'))
    g.declare_reviewed_scope(['tagmented_cdna','cas9','guide_pool','rnp_buffer','assemble_rnp','cut','proteinase','proteolyze','organic_extract_reagent','organic_extract','ethanol_recovery_reagent','ethanol_recovery','index_forward','index_reverse','kapa','index_pcr','size_selection_reagent','size_selection','dna_cleanup_reagent','dna_cleanup','q_pcr','other_libraries','pool','pooled_library'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_input_geometry','Exact NuGEN/Nextera input, full guide RNA transcription and historical index/primer/read geometry','Published Methods191/197/211 and directly cited originals'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();dna='ACGA'+GUIDE+'AGG'+'TCGA';original=Fragment(dna,frozenset({'left','right'}));cut=SelectedCas9Break(original,GUIDE,24);a,b=cut.products()
    c('selected cut lies3bases before PAM',cut.cut_position,21)
    c('both physical DNA products preserve complete sequence',a.dna+b.dna,original.dna)
    c('cut cannot create two original outer sites on either product',(original.paired_adapter_eligible,a.paired_adapter_eligible,b.paired_adapter_eligible),(True,False,False))
    try:SelectedCas9Break(Fragment(dna[:25]+'AA'+dna[27:],original.outer_sites),GUIDE,24);bad=False
    except ValueError:bad=True
    c('nonNGG is outside selected exact-match source relation',bad,True)
    source=Source(sources_root()/'dash__10.1186+s13059-016-0904-5'/'s13059-016-0904-5_supp_13059_2016_904_MOESM1_ESM.pdf.txt','published guide table')
    if have(c,source,label='original recognition-domain transcription'):c('selected guide sequence printed in source',GUIDE in source.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
