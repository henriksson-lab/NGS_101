"""NEB Q5/KLD control-example construct workflow, manual Version2.0_7/20."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Row, Scene, Segment, complement_segments, revcomp
from circular import CircularizationGeometry
from page import head, info
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import Drawing, closed_circle, role
from protocol_views import OligoEntry, OligoSection, WorkflowSection, assemble_page

PAGE_FILE = 'q5_sdm.html'
TITLE = 'Q5 site-directed mutagenesis — back-to-back PCR and KLD'
SOURCE_URL = 'https://www.neb.com/en/-/media/nebus/files/manuals/manuale0554.pdf?hash=B7F6605B20F0E0E661C478FF15479353&rev=53eb54a2178645ed8d2a1ac82c55afc8'
SOURCE = 'NEB E0554S instruction manual Version2.0_7/20, pp2–5 and7; control-example partial product only'
FORWARD = 'AACCCTGGCGttAcCcAACTTAATCG'
REVERSE = 'TTCCCAGTCACGACGTTGTAAAA'


def primer(name, sequence):
    return Construct((Segment(name, sequence),), name)


def intended_product():
    return Construct((Segment('source forward-primer written end', FORWARD, inferred=True),
        role('intervening retained control-plasmid interval; sequence/length unassigned', inferred=True),
        Segment('source reverse-primer written end', revcomp(REVERSE), inferred=True)),
        'INFERRED partial linear control amplicon')


def linear_product(*, phosphorylated=False):
    product = intended_product()
    bottom = complement_segments(tuple(product))
    sc = Scene()
    sc.strand('amplicon', tuple(product), label=product.name,
        mod5='phosphate' if phosphorylated else '')
    sc.anneal('complement', bottom, to='amplicon',
        pair=(bottom[-1].name, product.segments[0].name),
        unpaired=tuple(s.name for s in bottom if s.is_role_token()),
        mod5='phosphate' if phosphorylated else '')
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    template_parts = (role('source control plasmid; full sequence and methylation sites unassigned'),)
    parent = g.input('parent', 'source methylated control plasmid template', K.MOLECULAR,
        primary=True, state=closed_circle('source parental circular template', template_parts))
    f = g.input('forward', 'source control SDM forward primer', K.OLIGO,
        state=primer('control forward', FORWARD))
    r = g.input('reverse', 'source control SDM reverse primer', K.OLIGO,
        state=primer('control reverse', REVERSE))
    q5 = g.input('q5', 'Q5 Hot Start High-Fidelity master mix', K.CATALYST)
    amplicon = g.product('amplicon', 'INFERRED partial intended linear PCR product', K.MOLECULAR,
        state=linear_product(), chemical_state=intended_product())
    residual = g.product('residual_parent', 'parental-template carryover in PCR mixture', K.MOLECULAR, state=parent.state)
    g.react('inverse_pcr', 'Back-to-back primer exponential whole-plasmid PCR', S.PCR,
        [(parent,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(q5,R.CATALYST)], [amplicon,residual],
        geometry=PCRGeometry.unavailable(Scene(), reason='Manual Figure3 specifies outward/back-to-back primer geometry. The selected control full map and mismatch coordinates are unassigned; no exact circular primer landing is claimed. The source circle is drawn in its input state.'))
    kinase = g.input('kld_kinase', 'kinase component of the KLD mixture', K.CATALYST)
    ligase = g.input('kld_ligase', 'ligase component of the KLD mixture', K.CATALYST)
    dpni = g.input('kld_dpni', 'DpnI component of the KLD mixture', K.CATALYST)
    phosphorylated = g.product('phosphorylated', 'INFERRED intended phosphorylated PCR product', K.MOLECULAR,
        state=linear_product(phosphorylated=True), chemical_state=intended_product())
    g.react('phosphorylate', 'KLD: phosphorylate the PCR product ends', S.PHOSPHORYLATION,
        [(amplicon,R.TEMPLATE),(kinase,R.CATALYST)], [phosphorylated])
    circle = g.product('candidate', 'INFERRED intended circular control-product class', K.MOLECULAR,
        state=closed_circle('intended source control circle', tuple(intended_product()),
            closure_label='intended KLD last-to-first sealing'), chemical_state=intended_product())
    g.react('seal', 'KLD: intramolecular ligation of the PCR product', S.CIRCULARIZATION,
        [(phosphorylated,R.TEMPLATE),(ligase,R.CATALYST)], [circle],
        geometry=CircularizationGeometry.known('double-strand intramolecular ligation', statement='The manual establishes KLD intramolecular ligation of the phosphorylated double-stranded PCR product. Exact selected template and terminal base states remain unassigned; no cohesive/blunt sequence is invented.'))
    digested = g.product('digested_parent', 'digested parental-template class', K.MOLECULAR)
    g.react('digest_parent', 'KLD: digest parental-template carryover', S.DIGESTION,
        [(residual,R.TEMPLATE),(dpni,R.CATALYST)], [digested])
    g.mark_terminal(Terminal(digested.id,Disposition.UNRESOLVED,'Digested template remains in the KLD mixture; no separate physical cleanup or zero parental background is asserted.'))
    g.constructed(circle, 'Intended construct class before transformation and clone screening; partial source-qualified context, not a verified clone or full plasmid map.')
    g.declare_reviewed_scope(('parent','forward','reverse','q5','inverse_pcr','kld_kinase',
        'kld_ligase','kld_dpni','phosphorylate','seal','digest_parent','candidate'),SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(), (
        head(TITLE), '<div class="wrap">', '<h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_q5_sdm.html">Research notes</a></p>',
        info('<a href="'+SOURCE_URL.replace('&','&amp;')+'">NEB E0554 manual Version2.0_7/20</a>.'),
        info('INFERRED: the control-example product has the source primer-written ends and an unknown intervening interval. Full template, mismatch coordinates and end chemistry are unassigned. No generic plasmid length or exact primer landing is claimed.'),
        info('KLD is one mixture/incubation. The graph separates phosphorylation, ligation and parental-template digestion to explain the molecular relations; it does not prescribe separate incubations. Primers for Q5/KLD are back-to-back, unlike complementary overlapping QuikChange primers.'),
        OligoSection('Published control primers', (OligoEntry('forward, source mixed case retained',tuple(primer('forward',FORWARD))),
            OligoEntry('reverse',tuple(primer('reverse',REVERSE))))),
        WorkflowSection('Linear amplification and KLD circularization'),
        info('The endpoint precedes transformation, biological recovery and clone verification. Correct source-primer transcription does not establish the full PCR product, reaction yield or clone genotype.'),
        '</div>'))


def selftest():
    from checks import Check, Source, have
    from chemistry_validation import check_graph, Status
    from storage import sources_root
    check=Check()
    check('no contradictory claims in the source-scoped graph',
        [c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY], [])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'manuale0554.pdf.txt',
        'NEB E0554 official pinned manual; tools/get_sources.py --slugs q5-site-directed-mutagenesis --no-papers')
    if have(check,source,label='vendor control oligos and revision'):
        text=source.path.read_text()
        check('source printed revision is the reviewed manual','Version 2.0_7/20' in text,True)
        check('forward primer matches original mixed case',FORWARD in text,True)
        check('reverse primer matches original',REVERSE in text,True)
    check.report()
