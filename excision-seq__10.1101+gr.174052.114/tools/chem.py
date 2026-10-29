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
    c.report()

PAGE_FILE = 'excision-seq.html'
TITLE = 'Excision-seq — source-limited uracil pre/post-digestion branches'
NOTES = '01_excision-seq.html'
SOURCE = '<a href="https://doi.org/10.1101/gr.174052.114">Bryan et al. (2014)</a>; published PMC4158761, Results/Figure1 and Methods. Detailed Supplemental Methods remain uninspected.'
CAVEAT = 'Selected uracil branches only. UDG, abasic and gap classes are schematic source events without genomic coordinates, exact termini or conversion efficiency. Detailed Supplemental Methods, source-printed enzyme identity and complete library geometry remain unresolved. Broken DNA is not physically discarded simply because it cannot seed the represented PCR. UVDE/photolyase branches remain separate and unassembled.'
READOUT = 'Source names Illumina MiSeq/HiSeq2000. Predigestion read starts report the removed base immediately upstream; post-digestion coverage reports relative survival/amplification and is inversely related to modified-base content. Exact run/adapter/index geometry is unassembled.'
PREVIEWS = (('Predigestion uracil branch', 'pre_output'), ('Post-digestion uracil branch', 'post_output'))
@dataclass(frozen=True)
class ExcisionState:
    stage: str
    def __post_init__(self):
        if self.stage not in {'input','abasic','gap','post_mixture','pre_library','post_library'}:
            raise ValueError('unsupported Excision-seq state')
    def rows(self):
        if self.stage == 'input':
            return molecule('source uracil-containing DNA class', (role('genomic flank'),Segment('embedded deoxyuridine base symbol','U'),role('genomic flank'))).rows()
        if self.stage == 'abasic':
            return mixture(('abasic DNA', (role('genomic flank'),role('abasic DNA site; scaffold/end chemistry unassigned'),role('genomic flank'))),('excised nucleobase', (role('free uracil nucleobase; not a DNA polymer'),))).rows()
        if self.stage == 'gap':
            return mixture(('DNA cleavage population', (role('gap-containing DNA and possible released genomic fragments; complementary strand geometry unassigned'),)),('reaction products', (role('excised base/abasic-site products; retained until documented separation'),))).rows()
        if self.stage == 'post_mixture':
            return mixture(('PCR-competent class', (role('intact adapter-bearing DNA strands'),)),('interrupted class', (role('excision-interrupted adapter-bearing strands and small products; remain in tube'),))).rows()
        if self.stage == 'pre_library':
            return mixture(('source amplicon library roles',(role('unassembled source adapter end'),role('genomic fragment derived from excision-associated boundary'),role('opposite unassembled end'))),('unassigned carryover',(role('nonamplified DNA and reaction components; final workup not reconstructed'),))).rows()
        return mixture(('source amplicon library roles',(role('amplified intact original strand-derived genomic sequence'),)),('nonamplified physical molecules',(role('broken DNA, residual templates and reaction components; no physical disposal established'),))).rows()

def protocol_graph():
    g = ProtocolGraph(TITLE,authority='Published PMC4158761 Results/Figure1 and Methods; detailed recipes deferred')
    dna = g.input('dna','source genomic DNA with elevated uracil; selected pre/post comparison',K.MOLECULAR,state=ExcisionState('input'),primary=True)
    pre = g.product('pre','DNA allocated to source predigestion branch',K.MOLECULAR,state=dna.state)
    post = g.product('post','DNA allocated to source post-digestion branch',K.MOLECULAR,state=dna.state)
    g.react('branches','Represent source pre/post preparations separately',S.SPLIT,((dna,R.TEMPLATE),),(pre,post))
    udg = g.input('udg','source UDG; full original conditions deferred',K.CATALYST)
    ap = g.product('abasic','UDG-treated DNA and released nucleobases',K.MOLECULAR,state=ExcisionState('abasic'))
    g.react('base_excision','Remove uracil nucleobase with UDG',S.DIGESTION,((pre,R.TEMPLATE),(udg,R.CATALYST)),(ap,))
    endo = g.input('gap_enzyme','main-text “T4 endonuclease IV”; exact source identity requires Supplemental Methods',K.CATALYST)
    gaps = g.product('gaps','source single-base-gap/fragment population and excision products',K.MOLECULAR,state=ExcisionState('gap'))
    g.react('gap_creation','Source gap-producing enzyme treatment; preserve all physical products',S.DIGESTION,((ap,R.TEMPLATE),(endo,R.CATALYST)),(gaps,))
    pre_recipe = g.input('pre_recipe','source end repair/A tail/adapters/PCR; detailed Supplemental Methods unreviewed',K.CONDITION)
    pre_library = g.product('pre_output','source predigestion library plus unassigned physical carryover',K.LIBRARY,state=ExcisionState('pre_library'))
    g.react('pre_library_boundary','Predigestion library preparation boundary',S.MAGIC,((gaps,R.TEMPLATE),(pre_recipe,R.CONDITION)),(pre_library,),reason='Main source names the operations but full original conditions, separations, oligos and exact products are unreviewed. No workup or perfect excision selection invented.')
    shear = g.input('shear_condition','source mechanical shearing; detailed settings deferred',K.CONDITION)
    sheared = g.product('sheared','mechanically sheared post-digestion input',K.MOLECULAR,state=dna.state)
    g.react('shear','Shear before adapter preparation in post-digestion route',S.FRAGMENTATION,((post,R.TEMPLATE),(shear,R.CONDITION)),(sheared,))
    post_recipe = g.input('post_recipe','source end polishing and adapter ligation before excision; exact recipe deferred',K.CONDITION)
    adapted = g.product('post_adapted','source adapter-bearing fragments before enzyme selection',K.MOLECULAR,state=molecule('adapter-bearing roles',(role('source adapter'),role('genomic DNA with uracil patterns unspecified'),role('opposite source adapter'))))
    g.react('post_adapter_boundary','Prepare adapters before post-digestion selection',S.MAGIC,((sheared,R.TEMPLATE),(post_recipe,R.CONDITION)),(adapted,),reason='Detailed Supplemental Methods not read; no other library recipe or exact oligos imported.')
    post_enzymes = g.input('post_enzymes','source excision-repair enzyme treatment; full identities/conditions deferred',K.CATALYST)
    interrupted = g.product('post_interrupted','intact and excision-interrupted DNA coexist',K.MOLECULAR,state=ExcisionState('post_mixture'))
    g.react('post_digest','Interrupt modified-base-bearing strands; do not physically discard them',S.DIGESTION,((adapted,R.TEMPLATE),(post_enzymes,R.CATALYST)),(interrupted,))
    amplification = g.input('post_pcr_recipe','source PCR and library workup, exact primers/conditions unassembled',K.CONDITION)
    post_library = g.product('post_output','amplicons from intact strands plus unassigned physical carryover',K.LIBRARY,state=ExcisionState('post_library'))
    g.react('post_pcr_boundary','Amplify intact strands without treating digestion as physical purification',S.MAGIC,((interrupted,R.TEMPLATE),(amplification,R.CONDITION)),(post_library,),reason='PCR-competent intact strands are preferentially amplified; broken DNA remains physical material. Exact PCR and any actual workup await Supplemental Methods.')
    for product in (pre_library,post_library):
        g.mark_terminal(Terminal(product.id,Disposition.UNRESOLVED,'Source-positive Illumina library endpoint; exact polymers, reaction products and run geometry unassembled.'))
    g.declare_reviewed_scope(('dna','branches','udg','base_excision','gap_enzyme','gap_creation','pre_recipe','pre_library_boundary','shear_condition','shear','post_recipe','post_adapter_boundary','post_enzymes','post_digest','post_pcr_recipe','post_pcr_boundary'),g.authority)
    g.add_source_item(SourceItem('supplemental_methods','Detailed published Supplemental Methods and oligonucleotides','Main Methods114–118 explicitly assigns recipes/oligos to Supplemental Methods; current ZIP only article figures.'))
    return g
