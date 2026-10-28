"""Source-owned bounded model; no inherited historical library recipe."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Construct, Scene, Segment, panel, oligo, revcomp
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, SequencingPlan, SourceItem
from protocol_model import role, molecule
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
import seqprimers as sp

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('species_' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, x, state, *, library=False):
    condition = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident + '_dna', 'recovered DNA', K.LIBRARY if library else K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically removed fractions', K.WASTE)
    g.react(ident, source, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source-described physical separation; no complete-purity assumption.')
    return out

def render_page():
    from chemdraw import workflow_from_protocol_graph, workflow_panel
    from page import head, info
    import html
    g = protocol_graph()
    terminal = next(iter(g.sequencing))
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p class="research-notes"><a href="' + NOTES + '">Research notes</a></p>',
             info(SOURCE), info(SUMMARY), info('Private partial source model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source oligos</h2>', *oligos()])
    parts.extend(['<h2>Source-qualified preparation</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    for title, content in globals().get('TOPOLOGY_PANELS', ()):
        parts.extend(['<h2>' + html.escape(title) + '</h2>', content])
    sc = Scene(); sc.strand('final', list(FINAL_LIBRARY), label='source library roles')
    parts.extend([panel(sc.rows(), caption=FINAL_CAPTION), '<h2>Sequencing boundary</h2>',
                  '<p>' + html.escape(SEQUENCING_UNAVAILABLE) + '</p>',
                  '<p>' + html.escape(SEQUENCING_ENDING) + '</p>', '</main>'])
    return '\n'.join(parts)

def selftest():
    g = protocol_graph(); page = render_page()
    assert g.publication_problems(), 'Incomplete source boundary must block finished publication'
    assert 'data-protocol-complete=' not in page, 'Partial preview must not advertise completion'
    print('PASS source-qualified graph and rendering; symbolic kit boundaries remain explicit')

PAGE_FILE = 'dnase-seq-boyle-2008.html'
TITLE = 'DNase-seq — Boyle 2008 Solexa MmeI tag branch'
NOTES = '01_dnase-seq-boyle-2008.html'
SOURCE = '<a href="https://doi.org/10.1016/j.cell.2007.12.014">Boyle et al. (2008)</a>; associated accepted manuscript PMC2669738, Experimental Procedures DNase Sample Preparation and High-Throughput Sequencing.'
SUMMARY = 'Biotinylated linker I marks DNase ends; MmeI digestion, streptavidin enrichment, linker II ligation, linker PCR and PAGE yield the source Solexa library, read with the printed linker-I primer.'
CAVEAT = 'Source-limited preview: associated accepted manuscript used explicitly; the linked chemistry supplement remains inaccessible in this pass. DNase-end preparation, exact MmeI tag-end geometry, intervening wash/release conditions and full strand assembly are unresolved. Computed linker pairing is a sequence relation, not proof of genomic-end ligation chemistry. The 454 and array branches are not merged. No degenerate NN window is labelled a UMI or sample index.'
LINKER_I_TOP = 'ACAGGTTCAGAGTTCTACAGTCCGAC'
LINKER_I_BOTTOM = 'GTCGGACTGTAGAACTCTGAAC'
LINKER_II_TOP = 'TCGTATGCCGTCTTCTGCTTG'
LINKER_II_BOTTOM = 'CAAGCAGAAGACGGCATACGANN'
PCR_REVERSE = 'CAAGCAGAAGACGGCATACGA'
PCR_FORWARD = 'AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA'
SEQUENCING_PRIMER = 'CCACCGACAGGTTCAGAGTTCTACAGTCCGAC'
OLIGOS = (('Linker I top', LINKER_I_TOP, '5′ Bio-'), ('Linker I bottom', LINKER_I_BOTTOM, '5′ P-'), ('Linker II top', LINKER_II_TOP, '5′ P-'), ('Linker II bottom', LINKER_II_BOTTOM, '5′-'), ('PCR reverse', PCR_REVERSE, '5′-'), ('PCR forward', PCR_FORWARD, '5′-'), ('Sequencing primer', SEQUENCING_PRIMER, '5′-'))
def oligos():
    return [oligo(name, [Segment(name, seq)], five=mod, three='-3′') for name, seq, mod in OLIGOS]

def linker_i():
    # Fully source-written strands; the matching core and its placement are computed.
    core = revcomp(LINKER_I_BOTTOM)
    if not LINKER_I_TOP.endswith(core):
        raise ValueError('source linker I strands do not have the expected terminal core')
    overhang = LINKER_I_TOP[:-len(core)]
    sc = Scene()
    sc.strand('top', [Segment('linker I exposed end', overhang), Segment('linker I core', core)], label='Linker I top', mod5='biotin')
    sc.anneal('bottom', [Segment('linker I core complement', LINKER_I_BOTTOM)], to='top', pair=('linker I core complement','linker I core'), label='Linker I bottom', mod5='phosphate')
    sc.note('top', 'Computed antiparallel overlap from the two printed oligos; genomic ligation-end preparation not assigned.')
    return sc

def linker_ii():
    core = revcomp(LINKER_II_TOP)
    if not LINKER_II_BOTTOM.startswith(core):
        raise ValueError('source linker II strands do not have the expected paired core')
    sc = Scene()
    sc.strand('top', [Segment('linker II core', LINKER_II_TOP)], label='Linker II top', mod5='phosphate')
    sc.anneal('bottom', [Segment('linker II core complement', core), Segment('degenerate ligation end', LINKER_II_BOTTOM[len(core):])], to='top', pair=('linker II core complement','linker II core'), label='Linker II bottom')
    sc.note('bottom', 'NN retains source degeneracy; no identifier meaning asserted. Full MmeI substrate geometry remains unreviewed.')
    return sc
TOPOLOGY_PANELS = (('Printed linker oligos — computed overlap only', panel(linker_i().rows())), ('Linker II — computed overlap only', panel(linker_ii().rows())))
FINAL_LIBRARY = Construct((role('linker-I PCR-derived source end; exact assembly deferred'), role('DNase-adjacent genomic tag; source readout 20 bases'), role('opposite linker-II-derived source end')), name='Solexa DNase library roles')
FINAL_CAPTION = 'Exact oligo reagents are displayed above. Source-qualified library roles do not invent a fully assembled product through the unreviewed end-preparation/MmeI boundary.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'The source prints a linker-I sequencing primer and 20 bp Illumina tags. Complete final-library polymer and exact run-primer landing geometry remain unassembled; no source-unavailability inference about other run/index settings.'
SEQUENCING_ENDING = 'Read starts report DNase-adjacent genomic positions. Identical 5-prime starts on a strand are computationally collapsed after sequencing; 454 ~108bp tags are a separate source branch.'
READ_SPECS = (sp.read_spec('Read 1',20,'DNase-adjacent genomic tag','include'),)
def boundary_geometry(state, why):
    return AdapterLigationGeometry('unavailable', 'unavailable', (state,), why)

def protocol_graph():
    g = ProtocolGraph(TITLE, authority='Boyle2008 associated manuscript: Solexa branch, DNase Sample Preparation1303 and readout1309')
    nuclei = g.input('nuclei', 'intact primary human CD4+ T-cell nuclei', K.MOLECULAR, state=molecule('nuclear genomic DNA', (role('chromatin genomic DNA; sequence and accessible loci unspecified'),)), primary=True)
    dnase = g.input('dnase', 'DNase I; exact nuclei preparation/titration deferred to relevant supplement', K.CATALYST)
    cut = g.product('cut', 'DNase-digested genomic DNA; fragments still together', K.MOLECULAR, state=molecule('DNase-digested DNA', (role('genomic DNase end; exact source-end distribution unreviewed'),)))
    g.react('dnase_digest', 'Digest intact nuclei with DNase I', S.DIGESTION, ((nuclei,R.TEMPLATE),(dnase,R.CATALYST)), (cut,))
    linker1 = g.input('linker1', 'biotin/phosphorylated linker I oligos as printed', K.OLIGO, state=linker_i())
    ligation1 = g.input('ligation1_condition', 'source linker I ligation; intervening end preparation/conditions unresolved', K.CONDITION)
    tagged = g.product('tagged', 'source linker-I-tagged genomic DNA', K.MOLECULAR, state=molecule('tagged genomic role', (role('source biotin-linker I end; exact genomic junction unresolved'),role('genomic DNA'))))
    g.react('ligate_linker1', 'Ligate biotinylated linker I to DNase ends', S.LIGATION, ((cut,R.TEMPLATE),(linker1,R.STRUCTURAL),(ligation1,R.CONDITION)), (tagged,), geometry=boundary_geometry(linker_i(), 'Printed adapter polymers are known and paired by computation; source genomic-end preparation and compatible ligation geometry are not reviewed because relevant supplement remains unavailable.'))
    mme = g.input('mmei', 'MmeI; recognition site at source linker-I 3-prime end', K.CATALYST)
    cleaved = g.product('cleaved', 'MmeI cleavage mixture retains tagged DNA and untagged genomic products', K.MOLECULAR, state=mixture(('biotin-tagged DNase-adjacent DNA', (role('linker-I-bearing genomic tag; exact offsets/ends unreviewed'),)), ('other cleavage products', (role('untagged genomic DNA cleavage products still in tube'),))))
    g.react('mmei_digest', 'Digest with MmeI without premature physical removal', S.DIGESTION, ((tagged,R.TEMPLATE),(mme,R.CATALYST)), (cleaved,))
    streptavidin = g.input('streptavidin', 'Invitrogen streptavidin beads', K.SELECTION)
    retained = g.product('retained', 'bead-enriched DNase tag population', K.MOLECULAR, state=molecule('retained DNase tag', (role('biotin-linker-I-bearing genomic tag; noncovalently bead-bound'),)))
    wash = g.product('wash', 'physically separated affinity flowthrough/wash fractions', K.WASTE)
    g.react('capture', 'Enrich DNase ends on streptavidin beads', S.AFFINITY_CAPTURE, ((cleaved,R.TEMPLATE),(streptavidin,R.SELECTION)), (retained,wash))
    g.discard(wash, 'Physical affinity separation of unretained cleavage products; digestion alone did not remove these DNA molecules.')
    linker2 = g.input('linker2', 'phosphorylated/NN linker II oligos as printed', K.OLIGO, state=linker_ii())
    ligation2 = g.input('ligation2_condition', 'source linker II ligation; exact conditions/workup deferred', K.CONDITION)
    adapted = g.product('adapted', 'source linker-I/tag/linker-II product roles', K.MOLECULAR, state=molecule('adapter-bearing genomic tag', (role('source linker I'),role('DNase-adjacent genomic tag'),role('source linker II'))))
    g.react('ligate_linker2', 'Ligate linker II after streptavidin enrichment', S.LIGATION, ((retained,R.TEMPLATE),(linker2,R.STRUCTURAL),(ligation2,R.CONDITION)), (adapted,), geometry=boundary_geometry(linker_ii(), 'Printed adapter core and NN end retained; exact MmeI cleavage-end geometry and full junction unreviewed. Pairing is not an enzyme-product claim.'))
    forward = g.input('forward', 'printed linker-specific forward PCR primer', K.OLIGO, state=molecule('PCR forward primer',(Segment('PCR forward',PCR_FORWARD),)))
    reverse = g.input('reverse', 'printed linker-specific reverse PCR primer', K.OLIGO, state=molecule('PCR reverse primer',(Segment('PCR reverse',PCR_REVERSE),)))
    pcr = g.input('pcr_condition', 'source PCR; exact program/polymerase source supplement deferred', K.CONDITION)
    amplified = g.product('amplified', 'linker-PCR product; library roles only', K.MOLECULAR, state=FINAL_LIBRARY)
    g.react('pcr', 'Amplify using the printed linker-specific primers', S.PCR, ((adapted,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(pcr,R.CONDITION)), (amplified,), geometry=PCRGeometry.unavailable(adapted.state, reason='Primer polymers are verified; full source strand assembly through unreviewed end/MmeI boundary remains unassigned, so exact template sites are not asserted.'))
    library = cleanup(g,'page_cleanup','Source PAGE purification after PCR; precise source gel details deferred',amplified,FINAL_LIBRARY,library=True)
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':20},acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('nuclei','dnase','cut','dnase_digest','linker1','ligation1_condition','tagged','ligate_linker1','mmei','cleaved','mmei_digest','streptavidin','retained','wash','capture','linker2','ligation2_condition','adapted','ligate_linker2','forward','reverse','pcr_condition','amplified','pcr','page_cleanup_condition','page_cleanup_dna','page_cleanup_waste','page_cleanup'),g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Relevant original chemistry supplement remains unreviewed', 'PMC2669738 NIHMS106456-supplement-supplemental.pdf; generic retrieval gate, no chemical review'))
    return g

def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    source = Source(sources_root() / 'dnase-seq-boyle-2008__10.1016+j.cell.2007.12.014' / 'j.cell.2007.12.014_PMC2669738.html.txt', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC2669738/ ; get_sources.py --slugs dnase-seq-boyle-2008__10.1016+j.cell.2007.12.014')
    if have(c,source):
        text = source.path.read_text()
        for name, sequence, modification in OLIGOS:
            c(name + ' transcription from DNase Sample Preparation', sequence in text)
    protocol_graph(); render_page()
    c('Incomplete source boundary blocks publication', bool(protocol_graph().publication_problems()))
    c('Partial preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c.report()
