"""Circular TAR source variant, illustrated with prepared neutral reference maps."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Row, Scene, Segment
from assembly_geometry import CircularOverlapReference
from protocol_graph import ProtocolGraph, Disposition, Terminal, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import Drawing, closed_circle
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE = 'tar_circular.html'
TITLE = 'Circular TAR (1996) — illustrative capture reference'
SOURCE = 'Larionov et al.1996 doi10.1073/pnas.93.24.13925; published full text PMC19470, Introduction, Materials and Methods TAR Cloning Vectors, and Figure1 caption. Single recipient with repeats at opposite ends and yeast-mediated formation of circular YACs; not the earlier telomere/linear-YAC TAR variant.'
EXAMPLE = 'Invented short neutral prepared DNA and recipient with two explicit homology stand-ins. These are not human repeats, the published vector, a functioning YAC, or a validated capture design. The source circular topology and the intended toy base order are distinct claims.'
LEFT = 'ACGTCAGTACGATCGTGACTACGA'
RIGHT = 'TCGAGTACCTGACATGTCAGATCA'
TARGET = 'GATCTACGAGTC'
RECIPIENT = 'ATCGTAGCTCAG'


def references():
    target = Construct((Segment('left homology stand-in', LEFT, inferred=True),
                        Segment('neutral captured body', TARGET, inferred=True),
                        Segment('right homology stand-in', RIGHT, inferred=True)), 'invented prepared target interval')
    recipient = Construct((Segment('right homology stand-in', RIGHT, inferred=True),
                           Segment('neutral recipient body', RECIPIENT, inferred=True),
                           Segment('left homology stand-in', LEFT, inferred=True)), 'invented prepared recipient')
    return (('target', target, len(LEFT), len(RIGHT)),
            ('recipient', recipient, len(RIGHT), len(LEFT)))


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT,
                      representation='illustrative_example', example_description=EXAMPLE)
    refs = references()
    inputs = []
    for identifier, con, _, _ in refs:
        inputs.append(g.input(identifier, 'INFERRED — ' + con.name, K.MOLECULAR,
                              primary=True, state=Scene.duplex(tuple(con)), chemical_state=con,
                              inferred=True, source=EXAMPLE))
    yeast = g.input('yeast', 'source yeast recombination/recovery host; preparation unmodelled', K.BIOLOGICAL)
    geometry = CircularOverlapReference(refs, 'Declared toy overlap order only; these invented homology stand-ins do not establish source repeat compatibility, target specificity, DNA-end preparation or host chemistry.')
    con = geometry.reference('intended circular capture reference')
    final = g.product('final', 'INFERRED — intended circular capture reference', K.MOLECULAR,
                      state=closed_circle(con.name, tuple(con), closure_label='intended circle; host closure unmodelled'),
                      chemical_state=con, inferred=True, source=EXAMPLE)
    other = g.product('other', 'unresolved substrates, alternatives and cellular material', K.MOLECULAR,
                      state=Drawing((Row(chunks=[('Uncaptured DNA, surplus homology copies, residual recipient and alternative recombination products have unassigned physical fates. Target liberation and whole-genome material are outside this prepared-input example.', None, False)]),)))
    g.react('tar_capture', 'INFERRED — circular TAR host boundary', S.MAGIC,
            [(m, R.TEMPLATE) for m in inputs] + [(yeast, R.TEMPLATE)], [final, other],
            geometry=geometry,
            reason='Published circular-TAR topology is represented by prepared neutral maps. Intracellular end processing, recombination, DNA repair, covalent closure, cellular recovery and actual YAC maintenance remain unmodelled. The toy lacks functional YAC elements; exact genome cleavage and source vector maps are not supplied.')
    g.mark_terminal(Terminal(other.id, Disposition.UNRESOLVED, 'No physical purification, complete conversion, target specificity or yield is asserted.'))
    g.constructed(final, 'Conditional circular constructed-DNA reference; source-scale capture, functional YAC maintenance and sequencing are outside scope.')
    g.declare_reviewed_scope(('target', 'recipient', 'yeast', 'tar_capture', 'final'),
                            SOURCE + ' Exact bases, chosen distinct overlap strings and shortened prepared-input maps are invented; Other source branches, complete maps and linked-original review remain deferred.')
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(), (head(TITLE), '<div class="wrap"><h1>' + TITLE + '</h1>',
        '<p class="research-notes"><a href="01_tar_circular.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1073/pnas.93.24.13925">Larionov et al.1996</a>. INFERRED — short neutral prepared maps, homology stand-ins and intended reference order. This illustrates the single-recipient circular branch; it is not the paper’s vector, human repeats, source-scale capture or a functioning YAC.'),
        WorkflowSection('Prepared target and single recipient → yeast boundary → intended circular reference'), '</div>'))


def selftest():
    from checks import Check
    from chemistry_validation import check_graph, Status
    c = Check()
    claims = check_graph(protocol_graph())
    c('TAR host capture remains unresolved', any(x.target == 'tar_capture' and x.status is Status.UNRESOLVED for x in claims), True)
    c('reference comparison remains conditional', any(x.target == 'tar_capture' and x.status is Status.CONDITIONAL for x in claims), True)
    c('no contradictory declared claims', [x for x in claims if x.status is Status.CONTRADICTORY], [])
    c.report()
