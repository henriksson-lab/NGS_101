"""Protocol-owned source-qualified partial model; no completed-publication claim."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry

def state(label):
    sc=Scene();sc.strand('population',[Segment(label,'['+label+']',placeholder=True)]);return sc

def clean(g,key,current,procedure):
    reagent=g.input(key+'_reagent',procedure,K.SELECTION)
    kept=g.product(key+'_kept','Recovered '+procedure+'; exact molecular recovery/carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Physical separation fractions; no complete target specificity assumed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((current,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste));g.discard(waste,'Actual source physical separation only.');return kept
SOURCE='Dominissini et al.2016 DOI10.1038/nature16998; published SI Note2/3; associated accepted-manuscript PMC4842015 Methods144/160, final-VOR procedure applicability unverified.'
@dataclass(frozen=True)
class SelectedRNA:
    residues:tuple[str,...]
    def __post_init__(self):
        if not self.residues or any(x not in {'A','C','G','U','m1A','m6A'} for x in self.residues):raise ValueError('explicit supported RNA residues required')
    def dimroth(self,positions):
        positions=frozenset(positions)
        if any(type(i)is not int or not 0<=i<len(self.residues) or self.residues[i]!='m1A' for i in positions):raise ValueError('selected Dimroth path requires actual m1A sites')
        return replace(self,residues=tuple('m6A' if i in positions else x for i,x in enumerate(self.residues)))
    def copy_selected(self,choices=None):
        choices={} if choices is None else choices;comp={'A':'T','m6A':'T','C':'G','G':'C','U':'A'};bases=[]
        for i,x in enumerate(self.residues):
            if x=='m1A':
                if choices.get(i) not in {'A','C','G','T'}:raise ValueError('native m1A readthrough/mismatch must be caller selected, not universally predicted')
                bases.append(choices[i])
            else:bases.append(comp[x])
        return ''.join(reversed(bases))

COMMON_SCOPE=('rna','fragment_reagent','fragment','fragment_cleanup_reagent','fragment_cleanup','antibody_beads','ip','competitor','elute')
REARRANGEMENT_SCOPE=('alkali','dimroth','alkaline_cleanup_reagent','alkaline_cleanup')

def graph_for_branch(rearranged=True):
    g=ProtocolGraph('m1A-seq — published SI and associated AM-qualified route',authority=SOURCE,representation='published_instance')
    current=g.input('rna','Source total or polyA RNA alternative inputs, never an implicit physical pool',K.MOLECULAR,primary=True)
    reagent=g.input('fragment_reagent','Source RNA Fragmentation Reagents70°C15min to minimize unintended rearrangement',K.CONDITION)
    out=g.product('fragmented','RNA fragments and residual reagents retained before column',K.MOLECULAR)
    g.react('fragment','Associated AM fragmentation',S.FRAGMENTATION,((current,R.TEMPLATE),(reagent,R.CONDITION)),(out,));current=clean(g,'fragment_cleanup',out,'RNA Clean&Concentrator column')
    antibody=g.input('antibody_beads','Anti-m1A pre-coupled ProteinG Dynabeads,3h4°C',K.SELECTION)
    bound=g.product('bound','Selected bead-associated m1A RNA and possible background',K.MOLECULAR);flow=g.product('ip_flow','Separated unbound RNA/fractions; recovery unquantified',K.WASTE)
    g.react('ip','Source m1A immunoprecipitation',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(antibody,R.SELECTION)),(bound,flow));g.discard(flow,'Selected IP branch, not the separate stoichiometry microarray experiment.')
    competitor=g.input('competitor','Source free m1A mononucleoside molecular competitor; dose/packagingunassigned, notcovalentRNAconversion',K.MOLECULAR)
    out=g.product('eluted','Competitively eluted RNA plus possible m1A competitor/carryover',K.MOLECULAR)
    residual=g.product('residual_beads','Residual beads, possibly retained RNA and reagent background; partition unassembled',K.MOLECULAR)
    g.react('elute','Published SI Note2 competitive elution branch; PK/TRIzol alternative not merged',S.MAGIC,((bound,R.TEMPLATE),(competitor,R.STRUCTURAL)),(out,residual),reason='Exact molecular binding and separation geometry not assembled.');current=out
    g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual bead material disposal/recovery unassembled.'))
    if rearranged:
        alkali=g.input('alkali','Associated AM50mMNa2CO3/2mMEDTApH10.4,1h at source50°C or60°C alternatives; no single temperature chosen',K.CONDITION)
        out=g.product('rearranged','Selected partial m1A→m6A RNA rearrangement population; RNA residues retained, no assumption complete conversion',K.MOLECULAR)
        g.react('dimroth','Published SI Note2 alkaline rearrangement after IP before RT',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(alkali,R.CONDITION)),(out,));current=clean(g,'alkaline_cleanup',out,'RNA Clean&Concentrator after rearrangement')
    kit=g.input('kit','Actual historicalNEBNextUltra RNA library molecularsupplies; exactedition/oligos/strandprocessing/stockpackagingunreviewed',K.MOLECULAR)
    library=g.product('library','Source library population; full DNA polymers/read sites/indices unassembled',K.LIBRARY)
    g.react('library_prep','Associated AM source first-strand/library kit handoff',S.MAGIC,((current,R.TEMPLATE),(kit,R.STRUCTURAL)),(library,),reason='Named historical kit operation, full source edition and internal transitions unreviewed; no borrowed recipe.')
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1',),limitation='Associated AM HiSeq2500/TruSeqSRv3; cycles, full primer/index geometry and final-VOR applicability unverified.'))
    g.declare_reviewed_scope(COMMON_SCOPE+(REARRANGEMENT_SCOPE if rearranged else ()),SOURCE)
    g.add_source_item(SourceItem('pending_historical_kit','Historical NEBNext Ultra kit and complete adapter/read geometry','Associated AM Methods160; published SI Note2'))
    g.add_source_item(SourceItem('pending_vor','Final published Methods applicability','https://doi.org/10.1038/nature16998'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check
    c=Check();x=SelectedRNA(('C','m1A','G','U'));y=x.dimroth({1})
    c('selected rearrangement retains four RNA residues and N-methyl identity',y.residues,('C','m6A','G','U'))
    c('selected m6A copy is complementary firstDNA, not RNA-oriented report',y.copy_selected(),'ACTG')
    try:x.copy_selected();bad=False
    except ValueError:bad=True
    c('native m1A copying choice remains explicit',bad,True)
    c('partial rearrangement does not mutate untreated RNA',x.residues,('C','m1A','G','U'))
    c('both source comparison graphs retain publication boundaries',all(graph_for_branch(v).publication_problems() for v in (False,True)),True);c('Competitiveelution usesmolecular freecompetitor',graph_for_branch().materials['competitor'].kind,K.MOLECULAR);c.report()
if __name__=='__main__':selftest()
