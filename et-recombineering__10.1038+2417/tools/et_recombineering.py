"""Hypothetical neutral replacement reference; ET source identity only is verified."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Row, Scene, Segment
from protocol_graph import ProtocolGraph, Disposition, Terminal, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import Drawing, closed_circle
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE = 'et_recombineering.html'
TITLE = 'RecE/RecT ET — hypothetical neutral replacement reference'
SOURCE = 'Zhang et al.1998 doi10.1038/2417; defining published abstract PMID9771703. The abstract identifies RecE/RecT-mediated homologous recombination in E.coli, but does not establish the chosen substrate map or replacement branch.'
EXAMPLE = 'Invented short neutral DNA reference only, not a published ET experiment, functioning plasmid, or validated laboratory design. The substrate topology and replacement geometry are hypothetical pending the original methods. No lambda-Red components or later ET variants are inherited.'
LEFT = 'ACGTCAGTACGATCGT'
RIGHT = 'TCGAGTACCTGACATG'
RETAINED = 'GCTAGTCAGATCGACT'
OLD = 'GATCTACG'
NEW = 'ATCGTAGC'


def dna(name, seq):
    return Construct((Segment(name, seq, inferred=True),), name)


def references():
    """Build the intended replacement from named regions, never a host-path prediction."""
    left = Segment('left homology', LEFT, inferred=True)
    right = Segment('right homology', RIGHT, inferred=True)
    retained = Segment('retained neutral target body', RETAINED, inferred=True)
    old = Segment('displaced neutral interval', OLD, inferred=True)
    new = Segment('incoming neutral interval', NEW, inferred=True)
    target = Construct((left, old, right, retained), 'invented circular target reference')
    donor = Construct((left, new, right), 'invented linear donor reference')
    result = Construct((left, new, right, retained), 'intended replacement reference')
    return target, donor, result


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT,
                      representation='illustrative_example', example_description=EXAMPLE)
    target, donor, result = references()
    a = g.input('target', 'INFERRED — hypothetical circular target', K.MOLECULAR,
                primary=True, state=closed_circle(target.name, tuple(target)),
                chemical_state=target, inferred=True, source=EXAMPLE)
    b = g.input('donor', 'INFERRED — hypothetical linear donor', K.MOLECULAR,
                primary=True, state=Scene.duplex(tuple(donor)), chemical_state=donor,
                inferred=True, source=EXAMPLE)
    host = g.input('host', 'source RecE/RecT-dependent E.coli system; implementation unassigned', K.BIOLOGICAL)
    final = g.product('final', 'INFERRED — intended replacement reference', K.MOLECULAR,
                      state=closed_circle(result.name, tuple(result), closure_label='intended circle; host closure unmodelled'),
                      chemical_state=result, inferred=True, source=EXAMPLE)
    background = g.product('background', 'unresolved displaced DNA, parents and alternative products', K.MOLECULAR,
                           state=Drawing((Row(chunks=[('Displaced interval, donor homology copies, parental DNA and alternative products have unassigned physical fates; no detached intact fragment or complete conversion is asserted.', None, False)]),)))
    g.react('host_replacement', 'INFERRED — hypothetical ET replacement boundary', S.MAGIC,
            [(a, R.TEMPLATE), (b, R.STRUCTURAL), (host, R.TEMPLATE)], [final, background],
            reason='The original full methods have not been inspected. This two-sided replacement topology is an invented reference, not a verified branch of Zhang1998. RecE/RecT processing, strand exchange, displaced-strand chemistry, host repair, selection and recovery remain unmodelled; homology identity does not establish any of them.')
    g.mark_terminal(Terminal(background.id, Disposition.UNRESOLVED, 'Physical byproducts and competing species remain unassigned.'))
    g.constructed(final, 'Hypothetical reference endpoint for private draft review; no experimental ET construct or sequencing output is claimed.')
    g.declare_reviewed_scope(('target', 'donor', 'host', 'host_replacement', 'final'),
                            SOURCE + ' Only RecE/RecT system identity is verified; the full example is explicitly hypothetical.')
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(), (head(TITLE), '<div class="wrap"><h1>' + TITLE + '</h1>',
        '<p class="research-notes"><a href="01_et_recombineering.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/2417">Zhang et al.1998</a>. INFERRED — hypothetical substrate topology, homology strings and replacement map. Only the defining abstract was inspected; this draft does not establish that the original method uses this branch.'),
        WorkflowSection('Hypothetical target and donor → unmodelled ET boundary → intended replacement reference'), '</div>'))


def selftest():
    from checks import Check
    from chemistry_validation import check_graph, Status
    c = Check()
    claims = check_graph(protocol_graph())
    c('ET host replacement remains unresolved', any(x.target == 'host_replacement' and x.status is Status.UNRESOLVED for x in claims), True)
    c('no contradictory declared claims', [x for x in claims if x.status is Status.CONTRADICTORY], [])
    c.report()
