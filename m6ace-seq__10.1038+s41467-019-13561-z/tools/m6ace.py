"""m6ACE selected protection event; source-limited graph, not a complete library."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Segment, Scene, feature, revcomp
import illumina as il
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, Terminal, Disposition
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
    return Construct([Segment('RA5 RNA body',il.TRUSEQ_SMALL_RNA_RA5),Segment('constrained eight-base UMI','NNNNNNNB',placeholder=True,feature=feature('m6ace_umi','umi','degenerate'))])

def primer_scene():
    a=Adapter();sc=Scene();sc.strand('adapter',[Segment('adapter including terminal ddC',a.template)],mod5='rApp',mod3='ddC: no 3′OH')
    sc.anneal('primer',[Segment('RT primer',a.rt_primer())],to='adapter',pair=('RT primer','adapter including terminal ddC'))
    return sc

def protocol_graph():
    g=ProtocolGraph('m6ACE-seq — conditional exonuclease protection',authority=SOURCE,representation='published_instance')
    # Synthetic sequence tests are separate; source bulk substrate sequence is variable.
    state={'adapter3':Adapter(),'rna':'variable 5′p bulk RNA / 3′adapter with terminal ddC','selected_site':'crosslinked m6A; protected position retained','soluble_products':'none before digest'}
    prior=g.input('bound','Source bead-captured UV-crosslinked RNA after washes; upstream ligation and actual input split unmodelled',K.MOLECULAR,primary=True,chemical_state=state)
    x=g.input('xrn1','XRN1 1U, source LiCl buffer,37°C1h',K.CATALYST)
    dig=g.product('digested','Protected RNA plus soluble upstream digestion products; retained RNA starts at selected m6A',K.MOLECULAR,chemical_state={**state,'soluble_products':'upstream nucleotides retained until washing'})
    g.react('digest','Selected protected-substrate path; no complete digestion/protection efficiency asserted',S.DIGESTION,((prior,R.TEMPLATE),(x,R.CATALYST)),(dig,))
    washed=g.product('washed','Bead-associated protected RNA after documented XRN1 washes',K.MOLECULAR,chemical_state={**state,'soluble_products':'separated by washes'})
    waste=g.product('wash_waste','Physical XRN1 wash fractions including released nucleotide material',K.WASTE)
    g.react('wash','Source post-XRN1 bead washes',S.CLEANUP,((dig,R.TEMPLATE),),(washed,waste));g.discard(waste,'Physical wash separation')
    pk=g.input('proteinase','ProteinaseK 1mg/ml, SDS,50°C1.5h',K.CATALYST)
    el=g.product('eluate','Eluted RNA and proteolysis products; residual antibody/peptide chemistry unassigned',K.MOLECULAR,chemical_state={**state,'soluble_products':'separated by washes'})
    g.react('elute','ProteinaseK source elution; do not assert perfect removal of crosslink remnants',S.DEPROTEINIZATION,((washed,R.TEMPLATE),(pk,R.CATALYST)),(el,))
    clean=g.product('clean_rna','Ethanol-precipitated RNA; source proceeds to 5′UMI ligation and SSIII',K.MOLECULAR,chemical_state={**state,'soluble_products':'separated by washes'})
    ew=g.product('ethanol_waste','Physical precipitation reject',K.WASTE)
    g.react('precipitate','Source ethanol precipitation',S.CLEANUP,((el,R.TEMPLATE),),(clean,ew));g.discard(ew,'Physical precipitation reject')
    g.mark_terminal(Terminal(clean.id,Disposition.UNRESOLVED,'Upstream adapter ligation/IP/input split and downstream 5′ligation, RT/PCR geometry/index/read plan remain unassembled. NextSeqPE75 reported, not an acquired model output.'))
    g.declare_reviewed_scope(['bound','xrn1','digest','wash','proteinase','elute','precipitate'],SOURCE)
    return g

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
    c('blocked adapter nucleotide survives the modeled physical workup',protocol_graph().materials['clean_rna'].chemical_state['adapter3'].terminal,TerminalC())
    original=Source(sources_root()/'m6ace-seq__10.1038+s41467-019-13561-z'/'s41467-019-13561-z_supp_41467_2019_13561_MOESM4_ESM.xlsx.txt','Supplementary Data1')
    if have(c,original,label='source oligo workbook'):
        text=''.join(original.path.read_text().split())
        c('source adapter includes actual terminal ddC','/5rApp/'+a.body+'/3ddC/' in text,True)
        c('source RT primer matches computed full template',a.rt_primer() in text,True)
        c('source UMI preserves constrained B',(il.TRUSEQ_SMALL_RNA_RA5+'NNNNNNNB').replace('T','U') in text,True)
    c.report()
if __name__=='__main__':selftest()
