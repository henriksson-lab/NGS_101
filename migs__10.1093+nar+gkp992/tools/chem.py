"""Source-owned bounded model; no inherited historical library recipe."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Construct, Scene, Segment, panel, oligo, revcomp
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, SequencingPlan, SourceItem, Disposition
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
    parts.extend(['<h2>Prepare and physically wash the affinity reagent</h2>',
                  workflow_panel(workflow_from_protocol_graph(g, 'washed_conjugates')),
                  '<h2>Source-qualified genomic preparation</h2>',
                  workflow_panel(workflow_from_protocol_graph(g, terminal))])
    for title, content in globals().get('TOPOLOGY_PANELS', ()):
        parts.extend(['<h2>' + html.escape(title) + '</h2>', content])
    sc = Scene(); sc.strand('final', list(FINAL_LIBRARY), label='source library roles')
    parts.extend([panel(sc.rows(), caption=FINAL_CAPTION), '<h2>Sequencing boundary</h2>',
                  '<p>' + html.escape(SEQUENCING_UNAVAILABLE) + '</p>',
                  '<p>' + html.escape(SEQUENCING_ENDING) + '</p>', '</main>'])
    return '\n'.join(parts)

def selftest():
    from checks import Check
    g = protocol_graph(); page = render_page()
    c = Check()
    c('Incomplete source boundary blocks finished publication', bool(g.publication_problems()))
    c('Partial preview does not advertise completion', 'data-protocol-complete=' not in page)
    c('Genomic capture receives physically washed conjugates',
      [port.material for port in g.steps['capture'].inputs if port.role is R.SELECTION],
      ['washed_conjugates'])
    c('Conjugate wash has a physical discard terminal',
      g.terminals['conjugate_wash_waste'].disposition, Disposition.DISCARD)
    for material in ('captured', 'released', 'elution_cleanup_dna', 'library'):
        state = g.materials[material].state
        blocker_rows = [strand for strand in state.strands.values()
                        if any(seg.name == BLOCKER.name for seg in strand.segs)]
        c(material + ' retains possible blocker as a separate molecular species',
          len(blocker_rows) == 1 and len(blocker_rows[0].segs) == 1)
    c.report()

PAGE_FILE = 'migs.html'
TITLE = 'MiGS — MBD2 affinity isolation of methylated genomic fragments'
NOTES = '01_migs.html'
SOURCE = '<a href="https://doi.org/10.1093/nar/gkp992">Serre et al. (online 2009; issue 2010)</a>; published PMC2811030 Materials and Methods and gkp992_1.pdf.'
SUMMARY = 'Sonicated genomic DNA is captured by His-MBD2 domain on anti-His/protein-G beads, washed and released by SDS/Proteinase K, then prepared with the original Illumina ChIP-Seq kit for GAII 36 bp reads.'
CAVEAT = 'Partial source-owned capture model. Exact historical Illumina kit polymers, PCR/run primers and complete read architecture are unassembled. Supplemental bisulfite-validation primers are not MiGS library primers. Possible blocking-vector DNA carryover is retained separately through release and DNA recovery; neither proteolysis nor extraction/precipitation establishes its complete removal. Its recovery and eventual kit/read contribution are unassigned. The published elution buffer prints “20 ml Tris (pH 8)”; this wording is flagged rather than silently replaced with a concentration. Captured-fragment enrichment is not a base-resolution methylation call.'
BLOCKER = role('possible pCR4-TOPO blocking DNA; recovery and eventual kit fate unassigned')
FINAL_LIBRARY = Construct((role('unassembled original ChIP-Seq kit end'), role('MBD-enriched genomic sequence; ordinary DNA letters'), role('opposite source kit end')), name='MiGS GAII library roles')
FINAL_CAPTION = 'Representative source-positive genomic library roles; this panel does not establish removal of blocking DNA or exclude blocker-derived material. No other Illumina library recipe or validation primer is imported.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'GAII 36 bp end reads; exact historical kit adapter/run primer and index architecture unassigned.'
SEQUENCING_ENDING = 'Enrichment is interpreted over genomic fragments/windows. Identical mapped starts are computationally collapsed after sequencing, without a declared UMI. A read C is not itself a site-specific methylation call.'
READ_SPECS = (sp.read_spec('Read 1', 36, 'genomic sequence from MBD-enriched fragments', 'include'),)
def protocol_graph():
    g = ProtocolGraph(TITLE, authority='Serre et al. published PMC2811030 Materials and Methods38–52; supplement branch classification')
    dna = g.input('dna', 'purified genomic DNA after source Proteinase K/phenol/chloroform/ethanol extraction', K.MOLECULAR, state=molecule('genomic DNA', (role('genomic DNA with heterogeneous methyl-CpG pattern'),)), primary=True)
    sonic = g.input('sonication', 'source sonication in BW buffer to approximately 150–600 bp', K.CONDITION)
    fragmented = g.product('fragmented', 'sonicated genomic fragments and reaction components', K.MOLECULAR, state=mixture(('genomic DNA', (role('sonicated genomic DNA'),)), ('buffer components', (role('independent sonication-buffer components'),))))
    g.react('shear', 'Sonicate genomic DNA', S.FRAGMENTATION, ((dna, R.TEMPLATE), (sonic, R.CONDITION)), (fragmented,))
    clean = cleanup(g, 'fragment_cleanup', 'QIAquick cleanup, source excludes fragments smaller than 100 bp; 5 µg per capture', fragmented, dna.state)
    beads = g.input('protein_g_beads', '41.7 µl/ml protein G magnetic bead slurry', K.SELECTION)
    antibody = g.input('antibody', '10 µg/ml anti-His RGS antibody', K.SELECTION)
    mbd = g.input('mbd', '160 nM His-tagged recombinant MBD2 methyl-CpG-binding domain', K.SELECTION)
    blocker = g.input('blocker', '1 µg/ml self-ligated pCR4-TOPO blocking DNA', K.MOLECULAR, state=molecule('blocking vector DNA', (role('self-ligated pCR4-TOPO; not a genomic insert or index'),)))
    conjugation = g.input('conjugation_condition', 'BW buffer, overnight rotation at 4C', K.CONDITION)
    conjugates = g.product('conjugates', 'source antibody/proteinG/MBD bead preparation with blocker', K.MOLECULAR)
    g.react('prepare_capture', 'Prepare the source MBD-domain bead affinity reagent', S.BINDING, ((blocker, R.TEMPLATE), (beads, R.SELECTION), (antibody, R.SELECTION), (mbd, R.SELECTION), (conjugation, R.CONDITION)), (conjugates,), changes_molecule=False)
    conjugate_wash = g.input('conjugate_wash_condition', 'Wash conjugates with one volume cold 1x BW buffer; resuspend in one volume fresh 1x BW buffer', K.CONDITION)
    washed_conjugates = g.product('washed_conjugates', 'physically washed and resuspended MBD bead conjugates; residual blocker unquantified', K.SELECTION)
    conjugate_waste = g.product('conjugate_wash_waste', 'physically removed conjugate wash; reagent/blocker partition unquantified', K.WASTE)
    g.react('conjugate_wash', 'Wash and resuspend bead conjugates before genomic capture', S.CLEANUP,
            ((conjugates, R.TEMPLATE), (conjugate_wash, R.CONDITION)), (washed_conjugates, conjugate_waste))
    g.discard(conjugate_waste, 'Source physical conjugate wash before DNA addition; no complete blocker removal or quantitative partition asserted.')
    capture_condition = g.input('capture_condition', '5 ml source conjugates per 5 µg fragmented DNA; overnight rotation at 4C, five BW washes', K.CONDITION)
    captured = g.product('captured', 'bead-retained enriched DNA population', K.MOLECULAR, state=mixture(('retained genomic DNA', (role('MBD-affinity-enriched genomic DNA; methylation pattern unspecified'),)), ('possible blocking DNA carryover', (BLOCKER,)), ('source reagent species', (role('independent reagent-associated material; no complete-purity claim'),))))
    flowthrough = g.product('flowthrough', 'physically removed capture/wash fractions', K.WASTE)
    g.react('capture', 'Capture with MBD affinity and physically wash', S.AFFINITY_CAPTURE, ((clean, R.TEMPLATE), (washed_conjugates, R.SELECTION), (capture_condition, R.CONDITION)), (captured, flowthrough))
    g.discard(flowthrough, 'Source physical capture and five washes; unbound does not mean every fragment is unmethylated.')
    pk = g.input('elution', 'SDS/Proteinase K buffer, 60C 2h; printed “20 ml Tris (pH8)” remains unresolved', K.CATALYST)
    released = g.product('released', 'eluted DNA and digestion/reagent products', K.MOLECULAR, state=mixture(('DNA', (role('MBD-enriched genomic DNA'),)), ('possible blocking DNA carryover', (BLOCKER,)), ('reagent debris', (role('Proteinase K digestion products and SDS; present until workup'),))))
    g.react('elute', 'SDS/Proteinase K release without treating digestion as purification', S.DEPROTEINIZATION, ((captured, R.TEMPLATE), (pk, R.CATALYST)), (released,))
    isolated = cleanup(g, 'elution_cleanup', 'Phenol/chloroform extraction; isopropanol precipitation; resuspend in 45 µl water', released, mixture(('recovered genomic DNA', (role('MBD-enriched genomic DNA'),)), ('possible co-recovered blocking DNA', (BLOCKER,))))
    kit = g.input('kit', '10 ng DNA; original Illumina ChIP-Seq Sample Prep Kit, manufacturer protocol', K.CONDITION)
    library = g.product('library', 'source prepared GAII library; insert distribution near 120 bp', K.LIBRARY, state=mixture(('representative genomic library roles', tuple(FINAL_LIBRARY)), ('blocking DNA or kit-derived material; fate unresolved', (BLOCKER,))))
    g.react('library_prep', 'Original Illumina ChIP-Seq library preparation boundary', S.MAGIC, ((isolated, R.TEMPLATE), (kit, R.CONDITION)), (library,), reason='Named kit and sequencing endpoint source-positive. Exact historical kit processing, adapter members, primer polymers and amplification/run settings unreviewed; no family recipe or validation-primer substitution.')
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE, read_lengths={'Read 1':36}, acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('dna','sonication','fragmented','shear','fragment_cleanup_condition','fragment_cleanup_dna','fragment_cleanup_waste','fragment_cleanup','protein_g_beads','antibody','mbd','blocker','conjugation_condition','conjugates','prepare_capture','conjugate_wash_condition','washed_conjugates','conjugate_wash_waste','conjugate_wash','capture_condition','captured','flowthrough','capture','elution','released','elute','elution_cleanup_condition','elution_cleanup_dna','elution_cleanup_waste','elution_cleanup','kit','library','library_prep'), g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Original Illumina ChIP-Seq manufacturer protocol not reviewed', 'Published PMC2811030 Library preparation and sequencing explicitly refers to original manufacturer kit protocol'))
    return g
