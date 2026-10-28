"""Nm-seq terminal OED chemistry; historical library recipe remains unresolved."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition
SOURCE='Dai et al.2017 DOI10.1038/nmeth.4294; publisher Supplementary Protocol1 steps7–19; accepted manuscript mechanistic text71–77 (VOR applicability unreconciled)'

@dataclass(frozen=True)
class EndPath:
    rna: str
    methylated: frozenset[int]=frozenset()
    end3: str='OH'
    released: tuple[str,...]=()
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU'):raise ValueError('nonempty RNA required')
        if any(not 0<=p<len(self.rna) for p in self.methylated):raise ValueError('methylated site out of range')
        if self.end3 not in {'OH','dialdehyde','phosphate'}:raise ValueError('unknown terminal chemistry')
        if self.end3=='dialdehyde' and len(self.rna)-1 in self.methylated:raise ValueError('selected terminal Nm cannot have the vicinal-diol oxidation product')
    def oxidize(self):
        if self.end3!='OH':raise ValueError('vicinal-diol path requires dephosphorylated end')
        return self if len(self.rna)-1 in self.methylated else replace(self,end3='dialdehyde')
    def eliminate(self):
        if self.end3=='OH' and len(self.rna)-1 in self.methylated:return self
        if self.end3!='dialdehyde' or len(self.rna)==1:raise ValueError('selected elimination requires oxidized terminal ribose and surviving upstream RNA')
        return replace(self,rna=self.rna[:-1],end3='phosphate',released=self.released+(self.rna[-1],))
    def dephosphorylate(self):
        if self.end3=='dialdehyde':raise ValueError('eliminate before phosphatase')
        return replace(self,end3='OH')
    def cycle(self):return self.oxidize().eliminate().dephosphorylate()
    @property
    def ligatable3(self):return self.end3=='OH'
    def cleanup(self):return replace(self,released=())

def protocol_graph():
    g=ProtocolGraph('Nm-seq — bulk RNA terminal enrichment prefix',authority=SOURCE,representation='published_instance')
    prior=g.input('fragments','Source fragmented, Antarctic-phosphatase-treated bulk RNA; variable internal Nm locations',K.MOLECULAR,primary=True,chemical_state={'end3':'OH','Nm':'retained terminal Nm resists selected oxidation path'})
    scope=['fragments']
    for n in range(1,10):
        final=n==9
        condition=g.input(f'periodate_{n}','Periodate/lysine-HCl37°C30min; ethylene-glycol quench',K.CONDITION)
        oe=g.product(f'oe_{n}','Terminal Nm retains3′OH; susceptible selected N loses terminal residue, leaves3′phosphate; released products retained',K.MOLECULAR,chemical_state={'Nm_path':'3′OH, no shortening','susceptible_path':'shortened RNA3′p plus released residue products','yields':'unassigned'})
        g.react(f'oxidation_elimination_{n}','Final OE without dephosphorylation' if final else f'OED cycle{n}: source OE and quench',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(condition,R.CONDITION)),(oe,));scope += [f'periodate_{n}',f'oxidation_elimination_{n}'];prior=oe
        if not final:
            sap=g.input(f'rsap_{n}','rSAP2U37°C30min directly after quench',K.CATALYST)
            deph=g.product(f'deph_{n}','RNA3′OH population plus released residue products; recovery and conversion fractions unassigned',K.MOLECULAR)
            g.react(f'dephosphorylate_{n}','Dephosphorylate before next cycle',S.DEPHOSPHORYLATION,((prior,R.TEMPLATE),(sap,R.CATALYST)),(deph,));scope += [f'rsap_{n}',f'dephosphorylate_{n}'];prior=deph
        clean=g.product(f'clean_{n}','Column-recovered RNA; final3′OH/3′phosphate distinction retained' if final else 'Column-recovered RNA for next cycle',K.MOLECULAR)
        waste=g.product(f'waste_{n}','Physical RNA-column purification rejects and released residue products',K.WASTE)
        g.react(f'cleanup_{n}','Source RNA Clean and Concentrator purification',S.CLEANUP,((prior,R.TEMPLATE),),(clean,waste));g.discard(waste,'Physical column purification reject');scope += [f'cleanup_{n}'];prior=clean
    pnk=g.input('pnk','T4PNK3′phosphatase-minus plus ATP; publisher protocol37°C30min',K.CATALYST)
    ready=g.product('prepared','5′phosphorylated RNA retaining final3′OH/3′p distinction; conditional ligation eligibility, not perfect enrichment',K.MOLECULAR)
    g.react('phosphorylate5','Source5′phosphorylation; preserve blocked3′ends',S.PHOSPHORYLATION,((prior,R.TEMPLATE),(pnk,R.CATALYST)),(ready,));scope += ['pnk','phosphorylate5']
    g.mark_terminal(Terminal(ready.id,Disposition.UNRESOLVED,'Post-PNK cleanup and cited historical NEBNextE7330S library, adapters/PCR/read geometry unassembled. Publisher protocol30min and accepted-manuscript60min PNK timing are source-version comparison pending VOR reconciliation, not an established publication discrepancy.'))
    g.declare_reviewed_scope(scope,SOURCE);return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=EndPath('ACGU',frozenset({1}));q=p.cycle().cycle()
    c('two cycles expose internal Nm and conserve released sequence',(q.rna,q.released,q.end3),('AC',('U','G'),'OH'))
    c('terminal Nm halts selected progressive shortening',q.cycle(),q)
    unmodified=EndPath('ACGU').oxidize().eliminate()
    c('final OE without phosphatase leaves unligatable RNA3′p',(unmodified.rna,unmodified.end3,unmodified.ligatable3),('ACG','phosphate',False))
    c('physical cleanup removes released material only',unmodified.cleanup().rna,unmodified.rna)
    try:unmodified.oxidize();bad=False
    except ValueError:bad=True
    c('3′phosphate cannot enter next diol oxidation without repair',bad,True)
    original=Source(sources_root()/'nm-seq__10.1038+nmeth.4294'/'nmeth.4294_41592_2017_BFnmeth4294_MOESM3_ESM.pdf.txt','publisher Supplementary Protocol1')
    if have(c,original,label='original protocol checks'):
        text=' '.join(original.path.read_text().split())
        c('source repeats OED seven additional times','Repeat steps 7-12 seven more times.' in text,True)
        c('source explicitly excludes final dephosphorylation','not perform dephosphorylation (steps 10-11)' in text,True)
    c.report()
if __name__=='__main__':selftest()
