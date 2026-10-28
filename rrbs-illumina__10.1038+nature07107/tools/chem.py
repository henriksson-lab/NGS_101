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

PAGE_FILE = 'rrbs-illumina.html'
TITLE = 'RRBS — 2008 Illumina reduced representation bisulfite sequencing'
NOTES = '01_rrbs-illumina.html'
SOURCE = '<a href="https://doi.org/10.1038/nature07107">Meissner et al. (2008)</a>; associated accepted manuscript PMC2896277 Methods Summary; publisher Supplementary Figure S2.'
SUMMARY = 'MspI size selection precedes methylated Illumina adapter ligation, bisulfite conversion and PCR. Physical U is distinguished from the T subsequently copied and reported.'
CAVEAT = 'Partial model from the associated accepted manuscript and publisher supplement; complete VOR comparison, historical adapter polymers/end geometry, original kit workup, PCR and run primers remain unreviewed. Genomic base classes are schematic, not assigned genomic coordinates or quantitative conversion predictions. No other RRBS library recipe is imported.'
@dataclass(frozen=True)
class BisulfiteState:
    stage: str
    adapters: bool = False
    debris: bool = False
    def __post_init__(self):
        if self.stage not in {'input', 'converted', 'copied'}:
            raise ValueError('unknown RRBS molecular stage')
    def rows(self):
        sc = Scene()
        classes = [('unmodified C', 'C'), ('5mC', 'C')]
        for n, (label, base) in enumerate(classes):
            value = 'U' if self.stage == 'converted' and n == 0 else 'T' if self.stage == 'copied' and n == 0 else base
            parts = [role('genomic flank'), Segment('physical base ' + value, value), role('genomic flank')]
            if self.adapters:
                parts = [role('copied source adapter-derived end; chemistry unassembled' if self.stage == 'copied' else 'source methylated adapter; exact polymer unassembled')] + parts + [role('opposite source adapter-derived end')]
            sc.strand(str(n), parts, label=label + ' starting class')
            sc.note(str(n), 'Physical product: ' + ('5mC' if n == 1 and self.stage != 'copied' else value) + '; source idealized path, not conversion yield')
        if self.debris:
            sc.strand('carryover', [role('reaction components; not covalently joined to DNA')], label='same tube until cleanup')
        return sc.rows()
FINAL_LIBRARY = Construct((role('unassembled source Illumina PCR end'), role('RRBS genomic copy; U-derived T, retained-C-derived C'), role('opposite source PCR end')), name='RRBS library roles')
FINAL_CAPTION = 'Source-positive library components; exact original polymer assembly and run primer landing sites remain unassigned.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'Illumina 1G Genome Analyser, 36 bp end reads; complete original adapter/run-primer/index architecture not assembled in this pass.'
SEQUENCING_ENDING = 'Alignment against a converted reduced-representation reference permits strand-aware C/T interpretation. This graph does not establish perfect conversion or discrimination of every modified cytosine.'
READ_SPECS = (sp.read_spec('Read 1', 36, 'bisulfite-derived genomic copy', 'include'),)
def protocol_graph():
    g = ProtocolGraph(TITLE, authority='Meissner 2008 associated manuscript Methods Summary and publisher Supplementary Figure S2')
    dna = g.input('dna', '1–10 µg purified mouse genomic DNA', K.MOLECULAR, state=BisulfiteState('input'), primary=True)
    enzyme = g.input('mspi', '10–100 U MspI', K.CATALYST)
    cut = g.product('cut', 'MspI digested genomic DNA; all fragments still present', K.MOLECULAR, state=BisulfiteState('input', debris=True))
    g.react('digest', 'MspI digestion', S.DIGESTION, ((dna, R.TEMPLATE), (enzyme, R.CATALYST)), (cut,))
    selected = cleanup(g, 'select', 'Size selection: 40–120 bp and 120–220 bp fractions', cut, BisulfiteState('input'))
    repair = g.input('repair', 'source fill-in and 3-prime terminal-A extension reagents; members unreviewed', K.CATALYST)
    ended = g.product('ended', 'filled and 3-prime-A-extended DNA', K.MOLECULAR, state=molecule('A-extended genomic strand', (role('filled genomic DNA; opposite end also A-extended'), Segment('terminal A', 'A'))))
    g.react('end_processing', 'Fill in and add terminal 3-prime A', S.END_REPAIR, ((selected, R.TEMPLATE), (repair, R.CATALYST)), (ended,))
    ready = cleanup(g, 'end_cleanup', 'Phenol extraction and ethanol precipitation', ended, ended.state)
    adapter = g.input('adapter', 'pre-annealed Illumina adapter with 5mC instead of C; original polymer unresolved', K.OLIGO, state=molecule('source methylated adapter', (role('exact original methylated adapter; end geometry unreviewed'),)))
    kit = g.input('ligation_kit', 'source Illumina DNA preparation kit and protocol', K.CONDITION)
    ligated = g.product('ligated', 'methylated adapter-bearing genomic DNA', K.MOLECULAR, state=BisulfiteState('input', adapters=True, debris=True))
    geometry_scene = molecule('ligation boundary', (role('historical adapter/substrate geometry not reconstructed from inspected originals'),))
    geometry = AdapterLigationGeometry('unavailable', 'unavailable', (geometry_scene,), 'Source-positive pre-annealed methylated adapter ligation; exact original adapter end and kit revision remain unreviewed, not established absent.')
    g.react('ligate', 'Ligate methylated Illumina adapters', S.LIGATION, ((ready, R.TEMPLATE), (adapter, R.STRUCTURAL), (kit, R.CONDITION)), (ligated,), geometry=geometry)
    clean = cleanup(g, 'ligation_cleanup', 'QIAquick cleanup', ligated, BisulfiteState('input', adapters=True))
    bis = g.input('bisulfite', 'Qiagen EpiTect Bisulphite kit', K.CONDITION)
    converted = g.product('converted', 'chemical U-containing bisulfite products; original 5mC retained in source diagram', K.MOLECULAR, state=BisulfiteState('converted', adapters=True, debris=True))
    g.react('conversion', 'Bisulfite treatment; C to U precedes copying', S.CHEMICAL_CONVERSION, ((clean, R.TEMPLATE), (bis, R.CONDITION)), (converted,))
    kit_workup = g.input('kit_workup', 'EpiTect source kit workup; original procedural details deferred', K.CONDITION)
    pcr_ready = g.product('pcr_ready', 'source bisulfite-treated DNA supplied to PCR', K.MOLECULAR, state=BisulfiteState('converted', adapters=True, debris=True))
    g.react('workup_boundary', 'Original EpiTect workup boundary', S.MAGIC, ((converted, R.TEMPLATE), (kit_workup, R.CONDITION)), (pcr_ready,), reason='Named kit treatment is source-positive, but original kit workup not read. Carryover is not silently discarded and complete purity is not asserted.')
    forward = g.input('forward', 'original PCR forward primer; member unresolved', K.OLIGO, state=molecule('PCR forward primer', (role('historical PCR forward primer'),)))
    reverse = g.input('reverse', 'original PCR reverse primer; member unresolved', K.OLIGO, state=molecule('PCR reverse primer', (role('historical PCR reverse primer'),)))
    pcr = g.input('pcr_condition', 'preparative-scale PCR in source; exact program/reagents deferred', K.CONDITION)
    amplified = g.product('amplified', 'copied RRBS library', K.MOLECULAR, state=BisulfiteState('copied', adapters=True, debris=True))
    g.react('pcr', 'Preparative PCR copies U-derived sequence as T', S.PCR, ((pcr_ready, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (pcr, R.CONDITION)), (amplified,), geometry=PCRGeometry.unavailable(molecule('converted DNA', (role('source adapter-bearing bisulfite product'),)), reason='Original PCR primer members and full adapter/strand assembly unassigned.'))
    purified = cleanup(g, 'pcr_cleanup', 'QIAquick purification of PCR products', amplified, BisulfiteState('copied', adapters=True))
    library = cleanup(g, 'final', '4% NuSieve gel: 130–210 or 210–310 bp; QIAquick gel recovery', purified, FINAL_LIBRARY, library=True)
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE, read_lengths={'Read 1': 36}, acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('dna','mspi','cut','digest','select_condition','select_dna','select_waste','select','repair','ended','end_processing','end_cleanup_condition','end_cleanup_dna','end_cleanup_waste','end_cleanup','adapter','ligation_kit','ligated','ligate','ligation_cleanup_condition','ligation_cleanup_dna','ligation_cleanup_waste','ligation_cleanup','bisulfite','converted','conversion','kit_workup','pcr_ready','workup_boundary','forward','reverse','pcr_condition','amplified','pcr','pcr_cleanup_condition','pcr_cleanup_dna','pcr_cleanup_waste','pcr_cleanup','final_condition','final_dna','final_waste','final'), g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Published article and original library kit/workup not fully reviewed', 'Publisher VOR DOI10.1038/nature07107 and original kit instructions cited in Methods Summary'))
    return g
