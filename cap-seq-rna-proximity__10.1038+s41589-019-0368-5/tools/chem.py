"""Protocol-owned source-qualified partial graph; no finished-publication marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def clean(g,key,x,procedure):
    reagent=g.input(key+'_condition',procedure,K.CONDITION)
    out=g.product(key+'_kept','Recovered source preparation; exact molecular recovery/carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_separated','Actual separated source workup fractions; no ideal molecular partition assumed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((x,R.TEMPLATE),(reagent,R.CONDITION)),(out,waste));g.discard(waste,'Actual source physical workup only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True):
    ports=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION);ports.append((stock,R.CATALYST if catalyst else R.CONDITION))
    out=g.product(key+'_product',product or label+'; reaction remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,ports,(out,),reason=reason,changes_molecule=changes_molecule);return out

def unknown_adapter(label):
    sc=Scene();sc.strand('adapter',[Segment(label,'[source end geometry unreviewed]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),'Source operation known; exact original oligo/end/junction assembly unreviewed, no absence claim.')
SOURCE='Wang et al.2019 DOI10.1038/s41589-019-0368-5; published SI MOESM15 Fig1–7/13–18,especially Fig16 ordered preparation. Main Methods/whole library geometry unreviewed.'
@dataclass(frozen=True)
class SelectedRNA:
    sequence:str
    marks:tuple[str,...]=()
    def __post_init__(self):
        if not self.marks:object.__setattr__(self,'marks',('native',)*len(self.sequence))
        if not self.sequence or set(self.sequence)-set('ACGU') or len(self.marks)!=len(self.sequence) or any(m not in {'native','oxidized_class','PA_adduct','clicked_biotin'} or (m!='native' and b!='G') for b,m in zip(self.sequence,self.marks)):raise ValueError('selected G-origin productive classes only; no universal site or adduct-structure prediction')
    def advance(self,positions,from_mark,to_mark):
        if (from_mark,to_mark) not in {('native','oxidized_class'),('oxidized_class','PA_adduct'),('PA_adduct','clicked_biotin')}:raise ValueError('only explicitly selected source-qualified reaction classes permitted')
        positions=frozenset(positions)
        if any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='G' or self.marks[i]!=from_mark for i in positions):raise ValueError('actual supplied substrate chemical class required')
        return replace(self,marks=tuple(to_mark if i in positions else m for i,m in enumerate(self.marks)))

def protocol_graph():
    g=ProtocolGraph('CAP-seq RNA proximity — published SI Fig16 route',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Source miniSOG-expressing cells; compartment/probe/illumination conditions for exact sequencing samples unassembled',K.MOLECULAR,primary=True)
    x=op(g,'photo_label',cells,'Source miniSOG/blue-light/PA photo-oxidative RNA labeling',S.CHEMICAL_CONVERSION,reagent='source miniSOG-expressing cells, blue light and PA; not HaloTag/DBF or APEX',product='Selected PA-conjugated RNA plus oxidized/unmodified/background molecules; productive polymeric adduct identity unassembled')
    x=op(g,'extract',x,'Source Fig16 RNA extraction',S.MAGIC,reagent='original extraction reagents/partition unreviewed',reason='Preparation order established; exact extraction and remnant/physical separation fates unassembled.')
    x=op(g,'fragment',x,'Source Fig16 RNA fragmentation BEFORE click conjugation',S.FRAGMENTATION,reagent='source RNA-fragmentation conditions unreviewed',product='Fragmented modified/unmodified RNA and reaction remnants; no exact sizes/termini assigned')
    x=op(g,'click',x,'Source biotin-azide conjugation of extracted fragmented alkynyl RNA',S.CHEMICAL_CONVERSION,reagent='source click chemistry, exact sequencing-sample conditions unassembled',product='Selected clicked-biotin RNA adducts, residual PA/unmodified RNA and free reagent/products retained')
    beads=g.input('beads','Source streptavidin-coated beads',K.SELECTION)
    bound=g.product('bound','Source biotin-associated RNA fraction plus possible nonspecific/carryover species',K.MOLECULAR);other=g.product('capture_other','Source nonselected companion fraction; exact further workup/disposal unassembled',K.MOLECULAR)
    g.react('capture','Source Fig16 streptavidin enrichment',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'No uninspected washing/disposal recipe assumed.'))
    fa=g.input('formamide','Source95%formamide elution, time/temperature unreviewed',K.CONDITION)
    eluate=g.product('eluate','Source eluted RNA-containing preparation; exact adduct/bead/ligand partition unassembled',K.MOLECULAR);residual=g.product('residual','Companion bead-associated material after source elution; possible residual RNA/remnants',K.MOLECULAR)
    g.react('elute','Source Fig16 formamide elution',S.MAGIC,((bound,R.TEMPLATE),(fa,R.CONDITION)),(eluate,residual),reason='Reported elution and companion material; exact affinity dissociation, RNA conformation and adduct fate unassembled.');g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Source elution does not establish complete RNA recovery or residual disposal.'))
    primer=g.input('rt_primer','Source RT primer, exact body/member/end unreviewed',K.OLIGO);rt=g.input('rt_enzyme','Source RT enzyme/recipe unreviewed',K.CATALYST)
    cdna=g.product('cdna','Source RNA-derived cDNA plus original RNA/adduct/reagent remnants; copying outcomes unassigned',K.MOLECULAR)
    g.react('rt','Source Fig16 reverse transcription',S.REVERSE_TRANSCRIPTION,((eluate,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    prep=g.input('ngs_reagents','Source NGS library construction details blocked on unreviewed main; not inferred from qPCR/in-vitro validation',K.CONDITION)
    library=g.product('library','Source CAP-seq cDNA library; full ends/strands/primer/run geometry unassembled',K.LIBRARY)
    g.react('library_boundary','Opaque source NGS library continuation',S.MAGIC,((cdna,R.TEMPLATE),(prep,R.CONDITION)),(library,),reason='SI Fig16 establishes sequencing endpoint; actual downstream original assembly/cleanup/adapters remain unreviewed.')
    g.sequence(library,SequencingPlan(acquisition='unspecified',limitation='Source Fig16 NGS; exact acquired read roles/cycles/platform and priming sites unassembled.'))
    g.add_source_item(SourceItem('unreviewed_main','Defining main Methods and full original library/readout geometry','https://doi.org/10.1038/s41589-019-0368-5'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedRNA('ACGGU');o=x.advance({2},'native','oxidized_class');p=o.advance({2},'oxidized_class','PA_adduct');b=p.advance({2},'PA_adduct','clicked_biotin')
    c('productive selected chemistry preserves full RNA sequence',b.sequence,x.sequence)
    c('clicked adduct is chemically distinct from current free PA-alkyne class',b.marks,('native','native','clicked_biotin','native','native'))
    try:x.advance({2},'PA_adduct','clicked_biotin');bad=False
    except ValueError:bad=True
    c('click cannot bypass supplied alkynyl substrate',bad,True)
    c('oxidation alone does not assert all sites acquired PA',o.marks[2],'oxidized_class')
    g=protocol_graph();c('published fragment-before-click order retained',g.steps['click'].inputs[0].material,'fragment_product')
    c('source geometry boundary prevents finished publication',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
