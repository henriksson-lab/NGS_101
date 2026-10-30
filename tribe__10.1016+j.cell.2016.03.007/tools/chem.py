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
SOURCE='McMahon et al.2016 DOI10.1016/j.cell.2016.03.007; associatedAMPMC5027142 intro36–42,limits126–132,Methods140–150; final-VOR/detailedProcedures applicability unresolved.'
@dataclass(frozen=True)
class EditedRNA:
    original:str
    inosine_positions:frozenset[int]=frozenset()
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGU') or any(type(i)is not int or not 0<=i<len(self.original) or self.original[i]!='A' for i in self.inosine_positions):raise ValueError('selected RNA deamination requires actual A-origin positions; no binding/edit-efficiency prediction')
    @property
    def chemical_residues(self):return tuple('I' if i in self.inosine_positions else b for i,b in enumerate(self.original))
    def copy_selected(self):
        # Selected guanosine-like I readthrough class; actual distributions unassembled.
        return ''.join(reversed([{'A':'T','C':'G','G':'C','U':'A','I':'C'}[b] for b in self.chemical_residues]))

def graph_for_branch(branch='s2_rna'):
    if branch not in {'s2_rna','s2_nascent','sorted_neuron'}:raise ValueError('source-defined RNA preparation branches only')
    g=ProtocolGraph('TRIBE — '+branch,authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Source cells expressing selected RBP–ADAR catalytic-domain fusion; no universal RBP target specificity assigned',K.MOLECULAR,primary=True)
    current=op(g,'editing',cells,'Source in-vivo RNA A-to-I editing; S2 copper induction24h, neuronal expression context distinct',S.CHEMICAL_CONVERSION,reagent='source RBP–ADARcd expression, native RNA substrates',product='Selected inosine-bearing RNA, unedited/background RNA and original cell material; source inferred G calls are not physical RNA-G residues')
    if branch=='sorted_neuron':
        current=op(g,'manual_sort',current,'Source dissection/trituration and manual glass-micropipette sorting of fluorescent neurons',S.MAGIC,reagent='source Abruzzi2015 cell isolation procedure',reason='Source manual selection known; exact extraction/selection and molecular companion fates unreviewed. Fluorescence is not a sequencing barcode.')
    current=op(g,'rna_preparation',current,'Source '+('nascentRNA extraction' if branch=='s2_nascent' else 'RNA extraction for selected source population'),S.MAGIC,reagent='actual source extraction/predecessor recipe unreviewed',reason='Exact physical extraction/removal/recovery not supplied by a generic family recipe.')
    if branch=='s2_nascent':
        current=op(g,'rrna_removal',current,'Source-cited nascentRNA rRNA depletion',S.MAGIC,reagent='Pennington2013 source procedure, actual components unreviewed',reason='Named source operation; no guessed capture/digestion mechanism or discarded molecular fraction.')
        current=op(g,'polya_depletion',current,'Source two rounds of oligo-dT polyA depletion on nascentRNA preparation',S.MAGIC,reagent='source InvitrogenDynabeadsOligodT, historical protocol unreviewed',reason='Source depletion order known; exact bound/unbound/remnant partition and full nascent-library applicability unassembled.')
    kit=g.input('library_reagents','Source historical IlluminaTruSeq for S2 RNAseq' if branch!='sorted_neuron' else 'Source Abruzzi2015 low-input neuronal RNAseq preparation; unreviewed original',K.CONDITION)
    lib=g.product('library','Source RNA-derived sequencing library; exact primer/adapter/strand/end/workup geometry unassembled',K.LIBRARY)
    g.react('library_boundary','Opaque source library preparation continuation',S.MAGIC,((current,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='Actual source named kit/predecessor only; no present-day workflow, source-selected recipe or new UMI/cell identifier imposed.')
    g.sequence(lib,SequencingPlan(acquisition='unspecified',limitation='AssociatedAM establishes RNAseq; exact source acquisition/run/primer geometry and final-VOR applicability unresolved.'))
    g.add_source_item(SourceItem('gated_procedures','Actual detailed Procedures and original oligo tables','https://pmc.ncbi.nlm.nih.gov/articles/instance/5027142/bin/NIHMS767511-supplement-Procedures.pdf'))
    g.add_source_item(SourceItem('pending_predecessors','Historical library and direct extraction/depletion/sorting predecessors','AssociatedAMMethods142–148 cites Khodor2011,Pennington2013,Abruzzi2015 and IlluminaTruSeq; actual editions unreviewed.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check
    c=Check();x=EditedRNA('ACGA',frozenset({0}));first=x.copy_selected()
    c('physical RNA editing makes inosine, not guanosine',x.chemical_residues,('I','C','G','A'))
    c('selected firstDNA is complementary and RNA-oriented report contains G',(first,revcomp(first)),('TCGC','GCGA'))
    c('chemical marking preserves original RNA coordinates',len(x.chemical_residues),len(x.original))
    try:EditedRNA('ACGA',frozenset({1}));bad=False
    except ValueError:bad=True
    c('ADAR-selected A-to-I class cannot silently edit C as A',bad,True)
    c('all three source branches retain source/publication boundaries',all(graph_for_branch(b).publication_problems() for b in ('s2_rna','s2_nascent','sorted_neuron')),True)
    c('nascent polyA depletion is not a universal whole-S2 operation','polya_depletion' in graph_for_branch('s2_rna').steps,False);c.report()
if __name__=='__main__':selftest()
