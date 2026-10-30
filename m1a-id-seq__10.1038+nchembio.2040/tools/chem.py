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
SOURCE='Li et al.2016 DOI10.1038/nchembio.2040; published SI Fig3/5; main Methods unavailable. Published abstract scopes m1A-ID-seq; no complete library recipe claimed.'
@dataclass(frozen=True)
class SelectedRNA:
    residues:tuple[str,...]
    def __post_init__(self):
        if not self.residues or any(x not in {'A','C','G','U','m1A','m6A'} for x in self.residues):raise ValueError('explicit supported RNA chemical residues required')
    def demethylated(self,positions):
        positions=frozenset(positions)
        if any(type(i)is not int or not 0<=i<len(self.residues) or self.residues[i]!='m1A' for i in positions):raise ValueError('selected demethylation requires actual m1A sites')
        return replace(self,residues=tuple('A' if i in positions else x for i,x in enumerate(self.residues)))
    def selected_stop(self,index):
        if type(index)is not int or not 0<=index<len(self.residues) or self.residues[index]!='m1A':raise ValueError('selected m1A-mediated stop only, not a universal stall prediction')
        suffix=self.residues[index+1:]
        if any(x not in {'A','C','G','U'} for x in suffix):raise ValueError('other modified-site copying choices must be supplied, not predicted')
        return revcomp(''.join(suffix).replace('U','T'))

def protocol_graph():
    g=ProtocolGraph('m1A-ID-seq — published SI-qualified comparison',authority=SOURCE,representation='published_instance')
    rna=g.input('ip_rna','Published SI immunoprecipitated m1A-containing RNA; exact prior fragmentation/capture/elution recipe blocked on main Methods',K.MOLECULAR,primary=True)
    enzyme=g.input('demethylase','Source demethylation protein; exact enzyme identity/cofactors/temperature not established from selected SI alone',K.CATALYST)
    treated=g.product('treated','Selected demethylated RNA plus unconverted/background RNA; released methyl-derived products retained, identities unassembled',K.MOLECULAR)
    g.react('demethylate','Published SI Fig3 demethylation comparison, BSA omitted; exact recovery/order to RT unassembled',S.CHEMICAL_CONVERSION,((rna,R.TEMPLATE),(enzyme,R.CATALYST)),(treated,))
    primer=g.input('rt_primer','Source AMV priming reagent; exact oligo/body/ends unreviewed',K.OLIGO)
    rt=g.input('amv','Source AMV; SI Fig5b0.2mMdNTP/10mMMg optimized conditions',K.CATALYST)
    hybrid=g.product('hybrid','Original RNA retained with selected cDNA; RT-stop/readthrough distribution not predicted',K.MOLECULAR)
    g.react('rt','SI AMV comparison; copying does not alter physical RNA sequence',S.REVERSE_TRANSCRIPTION,((treated,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(hybrid,))
    g.mark_terminal(Terminal(hybrid.id,Disposition.UNRESOLVED,'Main Methods unavailable: exact intervening purification, sequencing-library assembly, adapters/PCR/run layout not read.'))
    g.declare_reviewed_scope(['ip_rna','demethylase','demethylate','rt_primer','amv','rt'],SOURCE)
    g.add_source_item(SourceItem('blocked_main','Defining main Methods and full library order','https://doi.org/10.1038/nchembio.2040'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedRNA(('C','m1A','G','U'))
    c('selected demethylation retains RNA positions and converts m1A to A',x.demethylated({1}).residues,('C','A','G','U'))
    c('native selected stop copies downstream suffix in firstDNA orientation',x.selected_stop(1),'AC')
    c('selected copy does not physically truncate RNA',x.residues,('C','m1A','G','U'))
    try:x.demethylated({0});bad=False
    except ValueError:bad=True
    c('demethylation cannot silently convert ordinary C',bad,True)
    c('source-blocked graph prevents finished publication',bool(protocol_graph().publication_problems()),True);c.report()
if __name__=='__main__':selftest()
