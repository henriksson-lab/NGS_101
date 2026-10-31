"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem, SequencingPlan
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('s' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, template, state):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.LIBRARY if template.kind is K.LIBRARY else K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source physical separation; no perfect recovery/depletion asserted.')
    return out

def end_boundary(label):
    sc = molecule('unassembled ligation-ready geometry', (role(label),))
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,), label + '; no source-absence claim.')

def render_page():
    import html
    g = protocol_graph()
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p><a href="' + NOTES + '">Research notes</a></p>', info(SOURCE),
             info('Private incomplete model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source ordering oligos</h2>', *oligos()])
    for label, terminal in PREVIEWS:
        parts.extend(['<h2>' + html.escape(label) + '</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    parts.extend([info(READOUT), '</main>'])
    return '\n'.join(parts)

def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Pending source/assembly boundary prevents finished publication', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())

    x=BrdUDNA('ATTC',frozenset({1}))
    c('IncorporatedBrdU remains a physical analogue',x.physical,('A','BrdU','T','C'))
    c('ReferenceTcall does not erasephysicalanalogue',x.reference_calls(),'ATTC')
    bad=False
    try:BrdUDNA('ATTC',frozenset({0}))
    except ValueError:bad=True
    c('Analogue cannot silentlylabelreferenceadenine',bad)
    c('Thymidinechase still permits physicalreplication',g.steps['chase'].changes_molecule,True)
    c('Actual BrdU/HU/pronase derive from grouped source supplies',all(g.producer_of(k).id=='label_components' for k in ('brdu','hu','pronase')),True)
    c('Thymidine chase contributes actual molecules',g.materials['chase_stock'].kind,K.MOLECULAR)
    c('Native cleanup preserves DNA library kind',g.materials['library'].kind,K.LIBRARY)
    c('Selected genomicbranch remains unsheared',all(s.kind is not S.FRAGMENTATION for s in g.steps.values()),True)
    c.report()
def pop(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other original DNA/RNA/protein or reaction remnants; exact partition and recovery unquantified'),)))
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,additional=()):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role))+tuple(additional),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))

PAGE_FILE='d-nascent.html'
TITLE='D-NAscent'
NOTES='01_d-nascent.html'
SOURCE='Associated manuscript EMS82269 Methods106/118/124; published SI MOESM1 selected genomic branch scope; 10.1038/s41592-019-0394-y.'
CAVEAT='Selected original yeast HU/BrdU-depletion genomic1DLSK109 route. Synthetic2D/1D2, barcoded training,8kb training shearing and BrdU-IP Illumina validation are separate. Historical end/adaptor assembly and final-VOR method applicability remain unreviewed.'
READOUT='Native MinIONMk1/R9.4 acquisition retains incorporatedBrdU chemical identity. HMM calls, segmentation and fork inference are computational; a calledT does not chemically remove the bromine or imply perfect incorporation/detection. No fixed native read length.'
PREVIEWS=(('SelectedHU native genomic branch','library'),)
from dataclasses import dataclass
@dataclass(frozen=True)
class BrdUDNA:
    sequence:str
    positions:frozenset[int]
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('explicit DNA sequence required')
        if not isinstance(self.positions,frozenset) or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='T' for i in self.positions):raise ValueError('BrdU replaces only actual reference thymidine positions')
    @property
    def physical(self):return tuple('BrdU' if i in self.positions else b for i,b in enumerate(self.sequence))
    def reference_calls(self):return self.sequence

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    x=g.input('yeast','source E3087 yeast population',K.MOLECULAR,state=pop('source yeastDNA with no caller-assigned BrdU positions'),primary=True)
    x=op(g,'arrest','Source alpha-factorG1 arrest',S.MAGIC,x,x.state,'source culture/arrest conditions')
    label_supplies=g.input('label_supplies','Actual HU/BrdU/pronase supply inventory; packaging and stock count unassigned',K.MOLECULAR)
    brdu=g.product('brdu','Actual40ug/ml BrdU molecular contribution',K.MOLECULAR)
    hu=g.product('hu','Actual200mM HU molecular contribution',K.MOLECULAR)
    pronase=g.product('pronase','Actual pronase-mediated release catalyst; dose unassembled',K.CATALYST)
    g.react('label_components','Expose actual source contributor classes; not a documented premix',S.GENERIC,((label_supplies,R.STRUCTURAL),),(brdu,hu,pronase),changes_molecule=False)
    x=op(g,'label','SelectedHU branch:40µg/mlBrdU plus200mMHU25min before pronaseS-release, then45minS-phase',S.MAGIC,x,pop('physical BrdU-bearing newly replicatedDNA plus unlabelledDNA/cellular components; incorporation positions/fractions unassigned'),'Source cell-replication/incubation environment; actual contributors supplied separately',additional=((brdu,R.STRUCTURAL),(hu,R.STRUCTURAL),(pronase,R.CATALYST)))
    x=op(g,'chase','Source400µg/mlthymidine chase;15minlater freshYPD+thymidine; replication continues, harvest times unresolved TableS1',S.MAGIC,x,x.state,'Actual400ug/ml thymidine-containing chase culture; nocodazoleG2/M timing unassembled',K.MOLECULAR,R.STRUCTURAL)
    x=op(g,'extract','Source genomicDNA extraction envelope:phenol/chloroform,RNaseA,proteinaseK,ethanol; exact internal order unassembled',S.MAGIC,x,pop('recovered BrdU-containingHMWgenomicDNA and possible carryover; extraction physical partitions unassembled'),'source genomic extraction conditions; no flow-cytometry SYTOX or BrdU-IP substitution')
    x=op(g,'repair','Source NEBE6050 end repair and genomic M6630FFPErepair; exact internal order/mark retention at repaired sites unreviewed',S.MAGIC,x,pop('source repaired HMWgenomicDNA; no added unsourced dA-tail; actual analogue retention at repaired sites unresolved'),'Actual E6050 and M6630 repair-supply inventory; packaging/component recipes unassembled, no premix asserted',K.MOLECULAR,R.STRUCTURAL)
    x=cleanup(g,'repair_cleanup','SelectedHU/depletion source0.4xAMPureXP cleanup',x,x.state)
    ds=mixture(('source1D adaptor',(role('originalLSK109 adaptor polymer and motor association; exact geometry unreviewed'),)))
    donor=g.input('adapter','source1DLSK109 adaptor',K.ADAPTER,state=ds);enzyme=g.input('ligase','source QuickT4DNA ligaseE6065',K.CATALYST)
    out=g.product('ligated','source1D adaptor-linked genomicDNA preparation',K.LIBRARY,state=pop('adaptor-linkedBrdU genomicDNA with associated motor; possible unligatedDNA/adaptor/reagents'))
    g.react('ligate','Source selected1D adaptor ligation',S.LIGATION,((x,R.TEMPLATE),(donor,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=AdapterLigationGeometry('unavailable','unavailable',(ds,),'Historical1DLSK109 genomic substrate/adaptor ends unreviewed; no source-absence claim.'))
    x=cleanup(g,'library','Source0.4xAMPureXP/proprietary long-molecule wash',out,out.state)
    native(g,x,'Source native MinIONMk1/R9.4 translocation, no sequencing primer or arbitrary nucleotide cycle count. Exact historical adapter/motor/end geometry unreviewed.')
    pending(g,'Original historical LSK109/E6050/M6630 kit documents, TableS1 sampling schedule and finalVOR applicability')
    return g
