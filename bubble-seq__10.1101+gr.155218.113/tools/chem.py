"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem
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
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
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
    c('Every DNA molecular output continues to an operation or explicit terminal',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c('Unselected source preparation fractions remain explicitly unresolved',g.terminals['isolation_other'].disposition is Disposition.UNRESOLVED)
    c('Bubble topology is not inherited as native topology of propagated clones',isinstance(g.materials['clones'].state,BubblePopulation) and g.materials['clones'].state.stage=='cloned')
    c.report()

def unassembled(label):
    return mixture(('source-derived DNA',(role(label),)),('possible carryover',(role('other original DNA/RNA, reaction fragments or protein carryover; quantities and exact partition unassigned'),)))

def operation(g,ident,label,kind,x,state,reagent,*,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    stock=g.input(ident+'_stock',reagent,reagent_kind)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular geometry/mixture relation not computed.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((x,R.TEMPLATE),(stock,reagent_role)),(out,),**kw)
    return out

def ligate(g,ident,label,x,adapter_label):
    adapter_scene=molecule(adapter_label,(role('source adapter; exact historical sequence/modifications and geometry unassembled'),))
    adapter=g.input(ident+'_adapter',adapter_label,K.ADAPTER,state=adapter_scene)
    enzyme=g.input(ident+'_enzyme','source T4DNA ligase',K.CATALYST)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=unassembled('source adapter-ligated product roles; no exact sealed junction claimed'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene,), 'Historical exact end/junction assembly remains unreviewed; no source-absence claim.')
    g.react(ident,label,S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=geometry)
    return out

def amplify(g,ident,label,x,primer_label):
    f=g.input(ident+'_forward',primer_label+' forward member; exact member unassigned',K.OLIGO,state=molecule('forward primer',(role('source primer identity unassembled'),)))
    r=g.input(ident+'_reverse',primer_label+' reverse member; exact member unassigned',K.OLIGO,state=molecule('reverse primer',(role('source primer identity unassembled'),)))
    enzyme=g.input(ident+'_enzyme','source Q5DNApolymerase',K.CATALYST)
    out=g.product(ident+'_product',label,K.LIBRARY,state=unassembled('source amplified library and residual original templates/cleaved species; no molecular deletion assigned to PCR'))
    g.react(ident,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(x.state,reason='Exact original adapters/primer/member/strand geometry unassembled.'))
    return out

PAGE_FILE='bubble-seq.html'
TITLE='Bubble-seq'
NOTES='01_bubble-seq.html'
SOURCE='<a href="https://doi.org/10.1101/gr.155218.113">Defining article</a>; Published PMC3814878 Introduction41–45, Results49 and Methods221–233; EuropePMC archive member inventory'
CAVEAT='Selected genomic replication-bubble route. Trapping/cloning precedents are explicitly source-named but unreviewed; nearby cDNA-library/RNA-seq instructions are a separate assay and are not imported. EuropePMC supplementary zip contains article figures only; actual publisher Supplemental Methods uninspected, web tool could not access its index. This is not evidence of absence.'
READOUT='Source pooled cloned genomic origin libraries sequenced on IlluminaGenomeAnalyzerII. Exact NGS construction, run and index geometry unresolved; read-depth origin calls and filters are computational.'
PREVIEWS=(("Source-qualified operation graph","library"),)

@dataclass(frozen=True)
class BubblePopulation:
    stage:str
    def __post_init__(self):
        if self.stage not in {'genomic','trapped','cloned','pooled','library'}:raise ValueError('unsupported bubble-seq state')
    def rows(self):
        sc=Scene()
        if self.stage in {'genomic','trapped'}:
            sc.strand('bubble',[role('EcoRI-bounded genomic replication-bubble restriction fragment class; branch/end topology unassembled')],label='replication-intermediate DNA class')
            sc.strand('other',[role('nonbubble genomic DNA and other original material; selected fraction not perfectly pure')],label='other physical material')
        else:
            sc.strand('clones',[role('pooled cloned origin-associated restriction-fragment DNA; host/vector propagation and exact cloned topology unassembled')],label='cloned genomic origin-associated DNA')
            sc.strand('other',[role('possible vector/host/background molecules; exact removal/partition unassigned')],label='other physical DNA')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    dna=g.input('genomic','source log-phaseGM06990 chromosomalDNA; ten RI isolations per duplicate library',K.MOLECULAR,state=BubblePopulation('genomic'),primary=True)
    isolation=g.input('ri_predecessors','source RI isolation/EcoRI fragment preparation and bubble trapping per named1991/2009/2011 predecessors; exact recipe unreviewed',K.CONDITION)
    trapped=g.product('trapped','source bubble-trapped genomic origin preparation',K.MOLECULAR,state=BubblePopulation('trapped'))
    other=g.product('isolation_other','physical RI/trap companion fractions; exact identity and disposition not reconstructed',K.MOLECULAR,state=unassembled('unselected companion DNA/other material; source-specific separation unreviewed'))
    g.react('isolate_trap','Source RI isolation and agarose bubble trapping; no invented enzyme/end or separation details',S.MAGIC,((dna,R.TEMPLATE),(isolation,R.CONDITION)),(trapped,other),reason='Original main defers preparation to named directly cited sources; source actual separation/fragment fates remain unassembled.')
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Do not assert molecular disposal without inspected preparation original.'))
    cloning=g.input('cloning','source cloning of origin DNA; actual host/vector/reaction details require chemistry originals',K.CONDITION)
    clones=g.product('clones','source cloned genomic origin library',K.MOLECULAR,state=BubblePopulation('cloned'))
    g.react('clone_boundary','Source cloning/propagation boundary; NGS rather than cDNA branch',S.MAGIC,((trapped,R.TEMPLATE),(cloning,R.CONDITION)),(clones,),reason='Main Results49 says after cloning; exact cloning chemistry/vector/propagation not inferred from other assays.')
    pooled=g.product('pooled','source pooled DNA from cloned origin libraries; separateB1/B2preparations retained as biological distinction',K.MOLECULAR,state=BubblePopulation('pooled'))
    pool=g.input('pool_condition','source pooled cloned libraryDNA, exact proportions unassigned',K.CONDITION)
    g.react('pool','Pool source cloned origin-library DNA for sequencing; no uniform clone recovery assumed',S.GENERIC,((clones,R.TEMPLATE),(pool,R.CONDITION)),(pooled,),changes_molecule=False,reason='Reported physical pooling, no new covalent product claimed.')
    kit=g.input('ngs_reagents','actual origin-library NGS construction pending SupplementalMethods; do not substitute adjacent cDNA GenomeAnalyzerIIx paragraph',K.CONDITION)
    library=g.product('library','source IlluminaGenomeAnalyzerII origin-DNA sequencing library',K.LIBRARY,state=BubblePopulation('library'))
    g.react('ngs_boundary','Unassembled genomic NGS preparation boundary',S.MAGIC,((pooled,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Primary source establishes NGS readout; actual origin-DNA NGS preparation is deferred, not imported from RNA-seq or vendor overview.')
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source GenomeAnalyzerII sequencing; molecular ends/adapters/run geometry unreviewed.'))
    g.add_source_item(SourceItem('supporting_methods','Actual publisher SupplementalMethods and directly named RI/trapping/cloning originals','https://genome.cshlp.org/content/suppl/2013/08/20/gr.155218.113.DC1 ; web tool could not access index, HTTP disposition unknown. EuropePMCzip contained only article figures.'))
    return g
