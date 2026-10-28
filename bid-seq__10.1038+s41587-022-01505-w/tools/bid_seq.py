"""BID-seq selected RNA-adduct/cDNA-deletion path, not a deletion-rate predictor."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Segment,feature,revcomp
import illumina as il
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition
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
    return Construct([Segment('3′ random barcode','NNNNN',placeholder=True,inferred=True,feature=feature('bid_random_3','umi','random',note='deduplication role inferred from source random-barcode processing')),
      Segment('inline barcode','ATCACG',feature=feature('bid_inline','inline_barcode','fixed')),
      Segment('3′adapter DNA',il.TRIM_SEEN_IN_READ1[:21]),Segment('non-nucleotide C3 block','[C3]',placeholder=True)])
def adapter5():
    return Construct([Segment('5′RNA body',il.TRUSEQ_SMALL_RNA_RA5),Segment('5′random barcode','NNNNN',placeholder=True,inferred=True,feature=feature('bid_random_5','umi','random',note='deduplication role inferred; distinct from opposite-end random barcode'))])

def graph_for_control(treated=True):
    g=ProtocolGraph('BID-seq — '+('treated' if treated else 'input'),authority=SOURCE,representation='published_instance')
    prior=g.input('ligated_rna','Both-adapter-ligated source RNA; 1.5µl input /8.5µl treated split is upstream',K.MOLECULAR,primary=True,chemical_state={'adapter3':adapter3(),'adapter5':adapter5(),'rna':'variable Ψ-containing RNA; RNA adduct is not a deletion'})
    scope=['ligated_rna']
    if treated:
        bs=g.input('bisulfite','2.4M Na2SO3/0.36M NaHSO3,70°C3h',K.CONDITION)
        changed=g.product('adduct_rna','Selected Ψ-BS RNA adduct path plus unassigned conversion/degradation fractions',K.MOLECULAR,chemical_state={'selected_site':'Ψ-BS','RNA_backbone':'site retained'})
        g.react('bisulfite_treat','Source BID-seq conditions; no universal reaction yield or C-deamination guarantee',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(bs,R.CONDITION)),(changed,));scope+=['bisulfite','bisulfite_treat']
        cleaned=g.product('desulfonated','Column-bound RNA after source desulfonation and washes',K.MOLECULAR,chemical_state={'selected_site':'Ψ-BS','RNA_backbone':'site retained'})
        waste=g.product('bs_waste','Physical bisulfite/column wash rejects',K.WASTE)
        g.react('desulfonate_clean','RNA Desulphonation Buffer1h on column; source washes/elution',S.CLEANUP,((changed,R.TEMPLATE),),(cleaned,waste));g.discard(waste,'Physical column rejects');prior=cleaned;scope+=['desulfonate_clean']
    primer=g.input('rt_primer','Printed SR RT primer',K.OLIGO,state=Construct([Segment('RT primer',RT_PRIMER)]))
    enzyme=g.input('ssiv','SuperScriptIV,50°C1h',K.CATALYST)
    hybrid=g.product('hybrid','RNA/cDNA hybrid with selected deletion and nondeletion paths; no fixed probability',K.MOLECULAR,chemical_state={'RNA_site':'Ψ-BS' if treated else 'Ψ','RNA_retained':True,'cDNA':'possible skipped complementary base only in treated selected path; background paths unmodelled'})
    g.react('rt','Source reverse transcription; RNA modification remains distinct from cDNA omission',S.REVERSE_TRANSCRIPTION,((prior,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST)),(hybrid,));scope+=['rt_primer','ssiv','rt']
    rnase=g.input('rnaseh','RNaseH,37°C20min',K.CATALYST)
    dig=g.product('rnase_digest','cDNA plus degraded RNA products retained in solution',K.MOLECULAR)
    g.react('rnase_h','Source RNA digestion; RNA products not yet physically discarded',S.DIGESTION,((hybrid,R.TEMPLATE),(rnase,R.CATALYST)),(dig,));scope+=['rnaseh','rnase_h']
    cdna=g.product('cdna','DNA-column-purified cDNA',K.MOLECULAR)
    waste=g.product('dna_waste','DNA cleanup rejects including RNA-digest material',K.WASTE)
    g.react('dna_cleanup','Source DNA Clean and Concentrator after70°C5min',S.CLEANUP,((dig,R.TEMPLATE),),(cdna,waste));g.discard(waste,'Physical DNA purification rejects');scope+=['dna_cleanup']
    g.mark_terminal(Terminal(cdna.id,Disposition.UNRESOLVED,'Exact historical NEBNext indexed PCR primer set/geometry and acquired NovaSeq6000SE100 read plan remain unassembled; upstream adapter ligations are also outside this graph.'))
    g.declare_reviewed_scope(scope,SOURCE);return g

def protocol_graph():return graph_for_control()
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
    c.report()
if __name__=='__main__':selftest()
