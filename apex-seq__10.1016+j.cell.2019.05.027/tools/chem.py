"""Protocol-owned source-qualified partial population graph, no finished marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def clean(g,key,x,procedure,retained='source-recovered preparation; exact molecular partition/carryover unquantified'):
    condition=g.input(key+'_condition',procedure,K.CONDITION)
    out=g.product(key+'_kept',retained,K.MOLECULAR)
    waste=g.product(key+'_separated','Source physically separated fractions; no perfect molecular recovery/removal claimed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(out,waste));g.discard(waste,'Actual source physical separation only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True):
    inputs=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION);inputs.append((stock,R.CATALYST if catalyst else R.CONDITION))
    out=g.product(key+'_product',product or label+'; reaction products/remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,inputs,(out,),reason=reason,changes_molecule=changes_molecule);return out
SOURCE='Fazal et al.2019 DOI10.1016/j.cell.2019.05.027; associated accepted-manuscript PMC6786773 Results93–99,Methods223/229/245/249. Final-VOR procedure applicability and gated published supplements unverified.'
@dataclass(frozen=True)
class SelectedDirectRNA:
    sequence:str
    guanine_adducts:frozenset[int]=frozenset()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU'):raise ValueError('explicit RNA required')
        if any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='G' for i in self.guanine_adducts):raise ValueError('this selected source G-adduct class requires actual G positions; broader site chemistry unmodelled')
    def label_selected(self,positions):
        return replace(self,guanine_adducts=self.guanine_adducts|frozenset(positions))
    def after_proteinase(self):
        return self # this selected substrate contains a direct RNA adduct, not a protein-borne label

def protocol_graph():
    g=ProtocolGraph('APEX-seq — AM-qualified polyA-selected route',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Source HEK293T cells expressing location-specific APEX2; construct/compartment branch explicit but not sequence-assembled',K.MOLECULAR,primary=True)
    loaded=op(g,'load_bp',cells,'BP500µM30min37°C',S.MAGIC,reagent='source biotinphenol/media',reason='Source exposure known; intracellular distribution/association unassembled.')
    tagged=op(g,'label',loaded,'APEX2/H2O2final1mM1min direct RNA-labeling population',S.CHEMICAL_CONVERSION,reagent='source H2O2, already expressed APEX2; not a new peroxidase stock',product='Selected covalent RNA/phenol-biotin adduct classes plus labeled proteins/unlabeled/background molecules; no perfect specificity')
    quenched=op(g,'quench',tagged,'Replace media with Trolox5mM/ascorbate10mM/azide10mMquench',S.QUENCHING,reagent='source quench solution')
    washed=clean(g,'cell_wash',quenched,'Source threeTrolox/ascorbateDPBS washes and harvest/pellet300g4min')
    rna=clean(g,'rna_extract',washed,'Source RNeasyPlusMini extraction with gDNAeliminator/RWTsubstitution; exact kit internals opaque','Source RNA preparation with direct RNA label and possible residual background; exact recovered polymers unassembled')
    beads=g.input('streptavidin','Source Pierce streptavidin magneticbeads2h4°C, then3B&Wwashes',K.SELECTION)
    bound=g.product('bound','Bead-associated biotinylatedRNA plus possible nonspecific/carryover material',K.MOLECULAR);flow=g.product('flowthrough','Source discarded supernatant/wash fractions; exact molecular recovery unquantified',K.WASTE)
    g.react('capture','Source RNA-biotin capture after RNA extraction',S.AFFINITY_CAPTURE,((rna,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.discard(flow,'Actual source magnet/supernatant discard and washes.')
    digested=op(g,'proteinase',bound,'ProteinaseK/sarcosine42°C1h then55°C1h',S.DEPROTEINIZATION,reagent='source ProteinaseK20mg/ml10µl plus digestionbuffer/Ribolock',catalyst=True,product='RNA-containing preparation with directRNA adducts, digested bead-associated proteins and other remnants; no automatic physical disposal or RNA-adduct reversal')
    purified=clean(g,'rna_cleanup',digested,'Source RNA Clean&Concentrator-5 purification')
    kit=g.input('truseq','Historical TruSeqstrandedmRNA kit including polyAselection; exact edition/oligos/internal transitions unreviewed',K.CONDITION)
    library=g.product('library','Source strandedpolyA-selected cDNA library; physical adapter/read/strand geometry unassembled',K.LIBRARY)
    g.react('library_prep','Source named TruSeqstrandedmRNA handoff; MITOtotalRNAalternative not merged',S.MAGIC,((purified,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Historical kit internals and molecular fates unassembled; no modern/predecessor recipe imported.')
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':75,'Read 2':75},limitation='AM sourceMiSeqQC thenHiSeq4000PE75; full library/index/primer geometry and final-VOR applicability unverified.'))
    g.add_source_item(SourceItem('pending_vor','Published final Methods and gated chemistry supplements','https://doi.org/10.1016/j.cell.2019.05.027'))
    g.add_source_item(SourceItem('pending_kit','Named historical extraction/TruSeq/RNAcleanup originals','AssociatedAMMethods229/245/249'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedDirectRNA('ACGGU').label_selected({2})
    c('direct selected RNA adduct preserves nucleotide positions',x.sequence,'ACGGU')
    c('direct RNA-bound label survives protein-only proteolysis relation',x.after_proteinase().guanine_adducts,frozenset({2}))
    try:SelectedDirectRNA('ACGGU').label_selected({1});bad=False
    except ValueError:bad=True
    c('selected G-adduct relation cannot silently relabel C as G',bad,True)
    c('unlabeled molecule need not be labeled merely because G-rich',SelectedDirectRNA('GGGG').guanine_adducts,frozenset())
    c('AM and kit boundaries prevent finished publication',bool(protocol_graph().publication_problems()),True);c.report()
if __name__=='__main__':selftest()
