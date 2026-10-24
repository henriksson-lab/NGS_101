"""SLIM insertion topology from Chiu et al.2004; no target-gene design."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S
from protocol_model import role, closed_circle, unsealed_circular_duplex
from protocol_views import WorkflowSection, assemble_page
from page import head, info

PAGE_FILE = 'slim.html'
TITLE = 'SLIM — four-primer inverse PCR and heteroduplex assembly'
SOURCE = 'Chiu et al.2004, doi:10.1093/nar/gnh172; Methods SLIM primer design, PCR and hybridization; Results Design of SLIM primers; original Figure1A; Table1 insertion tails'
INSERT = 'CACCATCACCATCACCAT'


def core():
    return role('retained plasmid body', inferred=True)


def tail(name='insertion tail'):
    return Segment(name, INSERT)


def duplex(sc, name, parts):
    bottom = complement_segments(parts)
    sc.strand(name, parts)
    sc.anneal(name + ' complement', bottom, to=name,
              pair=(bottom[0].name, parts[-1].name),
              unpaired=tuple(s.name for s in bottom if s.is_role_token()))


def products():
    sc = Scene()
    for name, parts in (
        ('Product1 FT + RS', (tail('left insertion'), core())),
        ('Product2 FS + RT', (core(), tail('right insertion'))),
        ('Product3 FT + RT', (tail('left insertion'), core(), tail('right insertion'))),
        ('Product4 FS + RS', (core(),))):
        duplex(sc, name, parts)
    return sc


def hybrids():
    """Figure1A productive hybrids; exact tail vs symbolic body are distinct."""
    sc = Scene()
    body = core()
    # Pair by the supplied body role only for layout, explicitly excluded from
    # molecular pairing checks. The unmatched tail is NOT paired to an N tract.
    top = (tail('5prime insertion'), body)
    bottom = complement_segments((body, tail('3prime insertion')))
    sc.strand('5prime-overhang hybrid top', top)
    sc.anneal('5prime-overhang hybrid bottom', bottom,
              to='5prime-overhang hybrid top', pair=(body.name + "'", body.name),
              unpaired=(body.name + "'",))
    top2 = (body, tail('3prime insertion'))
    bottom2 = complement_segments((tail('5prime insertion'), body))
    sc.strand('3prime-overhang hybrid top', top2)
    sc.anneal('3prime-overhang hybrid bottom', bottom2,
              to='3prime-overhang hybrid top', pair=(body.name + "'", body.name),
              unpaired=(body.name + "'",))
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    parent = g.input('parent', 'circular plasmid; gene body not assigned', K.MOLECULAR,
                     state=closed_circle('parent template', (core(),)), primary=True)
    primers = []
    for identifier, label, parts in (
        ('ft', 'FT: tailed forward', (tail(), role('FS binding core', inferred=True))),
        ('fs', 'FS: short forward', (role('FS binding core', inferred=True),)),
        ('rt', 'RT: tailed reverse', (Segment('reverse insertion tail', revcomp(INSERT)), role('RS binding core', inferred=True))),
        ('rs', 'RS: short reverse', (role('RS binding core', inferred=True),))):
        primers.append(g.input(identifier, label, K.OLIGO, state=Construct(parts, label)))
    pol = g.input('polymerases', 'Taq and Platinum Pfx DNA polymerases', K.CATALYST)
    dpni = g.input('dpni', 'DpnI', K.CATALYST)
    thermal = g.input('thermal', 'denaturation followed by two-stage cooling and hybridization', K.CONDITION)
    host = g.input('host', 'E.coli host repair and colony selection', K.SELECTION)
    pool = g.product('pcr_pool', 'INFERRED — four PCR product classes; target body schematic', K.MOLECULAR, state=products())
    g.react('inverse_pcr', 'Four-primer whole-plasmid inverse PCR', S.PCR,
            [(parent,R.TEMPLATE), *[(p,R.STRUCTURAL) for p in primers], (pol,R.CATALYST)], [pool],
            geometry=PCRGeometry.unavailable(Construct((core(),), 'circular template represented at a cut'),
                reason='This topology scope does not assign a target gene, binding-core sequence or full circular map. FT/FS share a core, RT/RS share the opposite core.'))
    cleaned = g.product('cleaned', 'INFERRED — product mixture after parental-template digestion', K.MOLECULAR, state=products())
    remnants = g.product('parent_remnants', 'digested parental DNA in the same tube', K.MOLECULAR, state=Construct((role('parental digest fragments'),),'parental digest fragments'))
    g.react('digest_parent', 'DpnI digest parental-template carryover', S.DIGESTION,
            [(pool,R.TEMPLATE),(dpni,R.CATALYST)], [cleaned,remnants])
    productive = g.product('productive', 'INFERRED — two productive heteroduplex classes', K.MOLECULAR, state=hybrids())
    other = g.product('other_hybrids', 'other reannealed classes; not expanded', K.MOLECULAR,
                      state=Construct((role('other heteroduplex population', inferred=True),),'other hybrids'))
    g.react('reanneal', 'Denature and reanneal the four products', S.THERMAL_TREATMENT,
            [(cleaned,R.TEMPLATE),(thermal,R.CONDITION)], [productive,other],
            reason='Figure1A: 16 possible duplex classes, two productive. This graph does not assign yields or equimolarity.')
    g.mark_terminal(Terminal(other.id, Disposition.UNRESOLVED,
        'The other hybrids remain in the reaction. Their subsequent molecular fates are not expanded.'))
    circle = g.product('nicked_circle', 'INFERRED — paired circular product; nicks remain', K.MOLECULAR,
        state=unsealed_circular_duplex((tail(),core()),
            statement='SLIM Figure1A: complementary exposed insertion tails hold the duplex in a noncovalent circle. The two nicks remain unsealed; their exact positions are not assigned in this schematic.'))
    g.react('pair_ends', 'Complementary exposed tails pair without in-vitro ligase', S.FOLDING,
            [(productive,R.TEMPLATE),(thermal,R.CONDITION)], [circle],
            reason='Original Figure1A and Results: noncovalent circles from productive hybrids. A limitation state preserves the unsealed duplex without inventing nick positions or a sealed strand; host repair is separate.')
    candidate = g.product('candidate', 'INFERRED — expected repaired insertion construct', K.MOLECULAR,
        state=closed_circle('expected host-repaired candidate', (tail(),core())))
    waste = g.product('unselected', 'material not carried into selected colonies', K.WASTE)
    g.react('host_recovery', 'Host repair and recovery of candidate colonies', S.MAGIC,
            [(circle,R.TEMPLATE),(remnants,R.TEMPLATE),(host,R.SELECTION)], [candidate,waste],
            reason='Host repair is not an in-vitro ligation. Colony recovery does not prove a particular genotype.')
    g.discard(waste,'Material not carried into recovered colonies; no quantitative fate or yield is assigned.')
    g.constructed(candidate,'Expected insertion topology; complete target sequence, clone genotype and experimental yield are outside scope.')
    g.declare_reviewed_scope(('parent','ft','fs','rt','rs','polymerases','dpni','host','thermal','inverse_pcr','digest_parent','reanneal','pair_ends','host_recovery','candidate'), SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(),(
        head(TITLE),'<div class="wrap">',f'<h1>{TITLE}</h1>',
        '<p class="research-notes"><a href="01_slim.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/nar/gnh172">Chiu et al., Nucleic Acids Research (2004)</a>; original Figure1A and Methods.'),
        info('Insertion branch. The published 18-base insertion tail is shown; target-specific primer cores and the retained plasmid are schematic inferred roles. No target gene, complete map or experimental clone is assigned. The other fourteen potential hybrids are not claimed to be physically removed.'),
        WorkflowSection('Four-primer assembly', paths=(('candidate',()),)),
        '<h2>Output</h2>',info('Constructed DNA. No sample index, UMI, NGS library or sequencing run is assigned. Deletion and substitution branches require separate target maps and are outside this insertion schematic.'),'</div>'))


def selftest():
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from storage import sources_root
    c=Check()
    c('no contradictory scoped molecular claims',[x for x in check_graph(protocol_graph()) if x.status is Status.CONTRADICTORY],[])
    s=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nar_gnh172_PMC535700.html.txt','Published article Table1 insertion primer tails')
    if have(c,s,label='published insertion-tail transcription'):
        text=s.path.read_text()
        c('forward tail in published tailed primer', INSERT+'GAGATTGACTTTACGGACGTA' in text)
        c('reverse tail in published tailed primer', revcomp(INSERT)+'AGGAGGGAAGTAATGTTTGTC' in text)
    c.report()
