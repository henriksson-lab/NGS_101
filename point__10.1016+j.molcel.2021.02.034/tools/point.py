"""Published POINT nascent-RNA isolation and source-specific library branch prefixes."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Sousa-Luís et al.2021 doi:10.1016/j.molcel.2021.02.034, VOR STAR Methods JATS487–503'
BRANCHES={'seq':('POINT-seq','NEBNext Ultra II Directional RNA E7760S'),
          '5':('POINT-5','SMARTer Stranded RNA-seq634836'),
          'nano':('POINT-nano','ONT direct cDNA SQK-DCS109')}

@dataclass(frozen=True)
class NascentRNA:
    carrier: str='RNA'
    artificial_poly_a: bool=False
    native_end_population: str='unassigned native RNA ends'
    def tail(self):
        if self.artificial_poly_a: raise ValueError('do not silently add the same artificial tail twice')
        return replace(self,artificial_poly_a=True)


def graph_for_branch(branch):
    if branch not in BRANCHES: raise ValueError('POINT branches are seq,5,nano')
    name,kit=BRANCHES[branch]
    g=ProtocolGraph(name+' — Pol-II-associated intact RNA',authority=SOURCE,representation='published_instance')
    sample=g.input('chromatin','nuclear/chromatin fraction from source cells;'+('4×10^7' if branch=='nano' else '10^7')+' cells',MaterialKind.BIOLOGICAL,primary=True)
    dnase=g.input('dnase','Turbo DNase in modified NUN/Empigen preparation',MaterialKind.CATALYST)
    digest=g.product('digested_chromatin','solubilized chromatin including RNA–PolII complexes and DNA digestion products',MaterialKind.MOLECULAR)
    g.react('digest_dna','Source modified NUN/Empigen chromatin preparation and DNase digestion',StepKind.DIGESTION,
            ((sample,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(digest,))
    soluble=g.product('soluble_fraction','soluble digested chromatin supernatant',MaterialKind.MOLECULAR)
    pellet=g.product('centrifuge_pellet','material rejected by source13,000rpm clarification',MaterialKind.WASTE)
    g.react('clarify','Centrifuge digested chromatin and retain supernatant',StepKind.CLEANUP,
            ((digest,InputRole.TEMPLATE),),(soluble,pellet));g.discard(pellet,'Actual clarification pellet')
    antibody=g.input('polii_beads','anti-PolII antibody-conjugated beads',MaterialKind.SELECTION)
    captured=g.product('captured_complex','PolII-associated intact RNA complexes',MaterialKind.MOLECULAR)
    washes=g.product('ip_waste','immunoprecipitation flow-through and washes',MaterialKind.WASTE)
    g.react('capture','Anti-PolII IP4°C1h; six NET-2E washes',StepKind.AFFINITY_CAPTURE,
            ((soluble,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(captured,washes));g.discard(washes,'Physical IP reject')
    state=NascentRNA()
    rna=g.product('purified_rna','Direct-zol-purified nascent RNA after further DNase treatment',MaterialKind.MOLECULAR,chemical_state=state)
    extract_waste=g.product('extraction_waste','extraction/column cleanup reject',MaterialKind.WASTE)
    g.react('extract','Source two Direct-zol purifications with intervening DNase treatment',StepKind.CLEANUP,
            ((captured,InputRole.TEMPLATE),),(rna,extract_waste));g.discard(extract_waste,'Physical extraction/column reject')
    endpoint=rna
    if branch=='seq':
        condition=g.input('fragmentation','historical E7760S fragmentation mixture;90°C10min',MaterialKind.CONDITION)
        endpoint=g.product('fragmented_rna','RNA-fragment mixture targeting150–200nt; all molecular fragments retained here',MaterialKind.MOLECULAR,
                           chemical_state=replace(state,native_end_population='mixture includes newly created fragment ends'))
        g.react('fragment','Fragment RNA before directional RNA library preparation',StepKind.FRAGMENTATION,
                ((rna,InputRole.TEMPLATE),(condition,InputRole.CONDITION)),(endpoint,))
    elif branch=='nano':
        pap=g.input('pap','E.coli poly(A) polymerase',MaterialKind.CATALYST)
        atp=g.input('atp','ATP required for poly(A) synthesis; source buffer recipe unreviewed',MaterialKind.MOLECULAR,inferred=True)
        tailed=g.product('poly_a_rna','RNA bearing in-vitro-added poly(A); tail length unassigned',MaterialKind.MOLECULAR,chemical_state=state.tail())
        g.react('tail','In-vitro polyadenylation before ONT direct-cDNA library prep',StepKind.TAILING,
                ((rna,InputRole.TEMPLATE),(pap,InputRole.CATALYST),(atp,InputRole.STRUCTURAL)),(tailed,))
        beads=g.input('size_beads','SPRISelect0.6×',MaterialKind.SELECTION)
        endpoint=g.product('selected_rna','source-selected polyadenylated RNA; removes RNA below500nt',MaterialKind.MOLECULAR,chemical_state=state.tail())
        small=g.product('small_rna_reject','RNA rejected by physical size selection',MaterialKind.WASTE)
        g.react('size_select','Source SPRISelect selection before SQK-DCS109',StepKind.SELECTION,
                ((tailed,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(endpoint,small));g.discard(small,'Physical size-selection reject')
    g.mark_terminal(Terminal(endpoint.id,Disposition.UNRESOLVED,
        kit+' downstream molecular inputs/oligos/revision and exact read geometry remain unreviewed. '
        +('POINT-5 uses random-N6/template switching without RNA fragmentation.' if branch=='5' else '')
        +' No other branch recipe is substituted.'))
    g.declare_reviewed_scope(('chromatin','dnase','digest_dna','clarify','polii_beads','capture','extract'),SOURCE)
    return g

def protocol_graph(): return graph_for_branch('seq')

def selftest():
    from checks import Check
    c=Check()
    c('all three source branches retain unresolved distinct endpoints',
      [sum(t.disposition is Disposition.UNRESOLVED for t in graph_for_branch(b).terminals.values()) for b in BRANCHES], [1]*3)
    c('POINT-nano polyadenylation remains RNA before direct-cDNA kit',
      graph_for_branch('nano').materials['selected_rna'].chemical_state.carrier,'RNA')
    c('POINT-5 prefix does not inherit fragmented-RNA branch', 'fragment' in graph_for_branch('5').steps,False)
    c.report()
if __name__=='__main__': selftest()
