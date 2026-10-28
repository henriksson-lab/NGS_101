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

PAGE_FILE = 'faire-seq.html'
TITLE = 'FAIRE-seq — formaldehyde/phenol enrichment with GAII readout'
NOTES = '01_faire-seq.html'
SOURCE = '<a href="https://doi.org/10.1016/j.ymeth.2009.03.003">Giresi and Lieb (2009)</a>; associated accepted manuscript PMC2710428, FAIRE and High-throughput sequencing Methods.'
SUMMARY = 'Recover nucleosome-depleted DNA preferentially in the aqueous phase of crosslinked chromatin, then prepare the source GAII A/T-overhang library.'
CAVEAT = 'Selected branch: human cultured cells and Illumina GAII; tissue, yeast, qPCR and microarray branches are not merged. Associated accepted manuscript used explicitly; one-page supplemental originals and VOR comparison deferred. Exact historical Illumina adapter/primer polymers remain unassembled. Fractionation is enrichment, not a binary molecular classifier of all regulatory sites.'
FINAL_LIBRARY = Construct((role('unassembled source Illumina PCR end'), role('FAIRE-enriched genomic sequence'), role('opposite original PCR end')), name='FAIRE GAII library roles')
FINAL_CAPTION = 'Source-positive GAII library roles; microarray LM-PCR linker sequences are not substituted for these ends.'
SEQ_PRIMERS = ()
SEQUENCING_UNAVAILABLE = 'GAII genomic end sequencing, 36 bp; original kit adapter/run primer/index members unassembled.'
SEQUENCING_ENDING = 'The paper computationally extends mapped 36 bp reads to 200 bp for density estimation. That operation does not synthesize DNA or acquire 200 sequenced bases.'
READ_SPECS = (sp.read_spec('Read 1', 36, 'FAIRE-enriched genomic end', 'include'),)
def chromatin(crosslinked=False, fragmented=False):
    sc = mixture(('genomic DNA', (role('fragmented genomic DNA' if fragmented else 'genomic DNA'),)), ('protein-associated genomic DNA', (role('protein-crosslinked genomic DNA' if crosslinked else 'chromatin-associated DNA'),)))
    sc.note('species_1', 'Source chromatin association; protein crosslinks are labelled, not base-resolved covalent structures.')
    return sc

def protocol_graph():
    g = ProtocolGraph(TITLE, authority='Giresi/Lieb 2009 associated manuscript: cultured human cell FAIRE and GAII sequencing Methods')
    cells = g.input('cells', 'human cultured-cell chromatin', K.MOLECULAR, state=chromatin(), primary=True)
    formal = g.input('formaldehyde', '1% formaldehyde, 5 min at room temperature', K.CONDITION)
    fixed = g.product('fixed', 'fixed chromatin mixture', K.MOLECULAR, state=chromatin(True))
    g.react('fix', 'Formaldehyde fixation of cultured human cells', S.CROSS_LINKING, ((cells, R.TEMPLATE), (formal, R.CONDITION)), (fixed,))
    gly = g.input('glycine', '125 mM glycine, 5 min', K.CONDITION)
    quenched = g.product('quenched', 'quenched fixed cells', K.MOLECULAR, state=chromatin(True))
    g.react('quench', 'Quench fixation', S.QUENCHING, ((fixed, R.TEMPLATE), (gly, R.CONDITION)), (quenched,))
    washed = cleanup(g, 'cell_wash', 'Source collection and two ice-cold PBS washes', quenched, chromatin(True))
    lysis = g.input('lysis_condition', 'source detergent/glass-bead lysis', K.CONDITION)
    lysed = g.product('lysed', 'chromatin lysate', K.MOLECULAR, state=chromatin(True))
    g.react('lysis', 'Lyse cells without inventing crosslink reversal', S.GENERIC, ((washed, R.TEMPLATE), (lysis, R.CONDITION)), (lysed,), changes_molecule=False, reason='Physical lysis releases chromatin; no covalent DNA transformation asserted.')
    sonic = g.input('sonication', 'Bioruptor high: 15 min, 30s on/off, 4C; source alternative kept separate', K.CONDITION)
    sheared = g.product('sheared', 'sonicated chromatin; bulk aliquot 100–1000 bp', K.MOLECULAR, state=chromatin(True, True))
    g.react('shear', 'Sonicate fixed chromatin', S.FRAGMENTATION, ((lysed, R.TEMPLATE), (sonic, R.CONDITION)), (sheared,))
    cleared = cleanup(g, 'clear', '15000g, 5 min, 4C; transfer supernatant from cellular debris', sheared, sheared.state)
    solvent = g.input('solvent', 'phenol/chloroform; repeat extraction then chloroform', K.SELECTION)
    aqueous = g.product('aqueous', 'aqueous FAIRE-enriched DNA fraction', K.MOLECULAR, state=mixture(('aqueous DNA', (role('genomic DNA preferentially depleted of nucleosome crosslinks'),)), ('RNA carryover', (role('independent RNA species'),))))
    nonaqueous = g.product('nonaqueous', 'interphase and other non-recovered fractions', K.WASTE)
    g.react('extract', 'Recover aqueous DNA, separating crosslinked material preferentially to interphase', S.EXTRACTION, ((cleared, R.TEMPLATE), (solvent, R.SELECTION)), (aqueous, nonaqueous))
    g.discard(nonaqueous, 'Physical phase separation; no complete-purity or perfect chromatin-class partition assumption.')
    recovered = cleanup(g, 'precipitate', 'Sodium acetate/glycogen/ethanol precipitation; 70% ethanol wash', aqueous, aqueous.state)
    rnase = g.input('rnase', 'RNase A, 1h at 37C', K.CATALYST)
    digested = g.product('rnase_treated', 'DNA and RNA digestion products remain together', K.MOLECULAR, state=mixture(('DNA', (role('FAIRE genomic DNA'),)), ('RNA debris', (role('RNA digestion products; physically retained until next cleanup'),))))
    g.react('rna_digest', 'Digest RNA without prematurely discarding its products', S.DIGESTION, ((recovered, R.TEMPLATE), (rnase, R.CATALYST)), (digested,))
    dna = cleanup(g, 'dna_cleanup', 'Source spin-column cleanup capable of recovering 75–200 bp DNA', digested, molecule('FAIRE DNA', (role('FAIRE-enriched genomic fragment'),)))
    repair = g.input('end_repair_reagent', 'Epicentre ER0720; blunt fragment ends', K.CATALYST)
    blunt = g.product('blunt', 'blunt-ended FAIRE DNA', K.MOLECULAR, state=molecule('blunt DNA role', (role('blunt-ended genomic fragment'),)))
    g.react('end_repair', 'Blunt FAIRE DNA ends', S.END_REPAIR, ((dna, R.TEMPLATE), (repair, R.CATALYST)), (blunt,))
    tail = g.input('a_tail_reagent', 'Epicentre KL06041K; add A overhang', K.CATALYST)
    tailed = g.product('tailed', '3-prime A-tailed DNA', K.MOLECULAR, state=molecule('A-tailed genomic strand', (role('source genomic body; opposite strand also A-tailed'), Segment('terminal A', 'A'))))
    g.react('a_tail', 'Add terminal A overhang', S.TAILING, ((blunt, R.TEMPLATE), (tail, R.CATALYST)), (tailed,))
    geometry = reported_adapter_end('Illumina #1000521 adapter', kind='cohesive', polarity='3-prime', terminal_sequence='T', statement='Source explicitly states double-stranded adapter with T overhang. Its full original body and complete duplex geometry are not assembled; source A overhang shown separately.', substrate_scenes=(tailed.state,))
    adapter = g.input('adapter', 'Illumina #1000521 double-stranded T-overhang adapters', K.OLIGO, state=geometry.scenes[0])
    ligase = g.input('ligation_condition', 'source GAII adapter ligation; exact historical formulation unreviewed', K.CONDITION)
    ligated = g.product('ligated', 'adapter-bearing FAIRE DNA', K.MOLECULAR, state=molecule('ligated library roles', (role('source adapter; body unassembled'), role('FAIRE genomic fragment'), role('opposite adapter'))))
    g.react('ligate', 'Ligate source T-overhang adapters to A-tailed DNA', S.LIGATION, ((tailed, R.TEMPLATE), (adapter, R.STRUCTURAL), (ligase, R.CONDITION)), (ligated,), geometry=geometry)
    selected = cleanup(g, 'gel_select', '2% agarose gel, excise 125–275 bp ligation products', ligated, ligated.state)
    forward = g.input('forward', 'Illumina #1000537 primer; original polymer unassembled', K.OLIGO, state=molecule('forward primer', (role('original #1000537 polymer'),)))
    reverse = g.input('reverse', 'Illumina #1000538 primer; original polymer unassembled', K.OLIGO, state=molecule('reverse primer', (role('original #1000538 polymer'),)))
    pfu = g.input('pfu', 'PfuUltraII, Stratagene #600670', K.CATALYST)
    library = g.product('library', 'PCR-amplified GAII library', K.LIBRARY, state=FINAL_LIBRARY)
    g.react('pcr', 'Amplify selected ligation products with source adapter primers', S.PCR, ((selected, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (pfu, R.CATALYST)), (library,), geometry=PCRGeometry.unavailable(selected.state, reason='Original #1000537/#1000538 polymers and complete adapter assembly unreviewed; microarray LM-PCR oligos excluded.'))
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE, read_lengths={'Read 1':36}, acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(('cells','formaldehyde','fixed','fix','glycine','quenched','quench','cell_wash_condition','cell_wash_dna','cell_wash_waste','cell_wash','lysis_condition','lysed','lysis','sonication','sheared','shear','clear_condition','clear_dna','clear_waste','clear','solvent','aqueous','nonaqueous','extract','precipitate_condition','precipitate_dna','precipitate_waste','precipitate','rnase','rnase_treated','rna_digest','dna_cleanup_condition','dna_cleanup_dna','dna_cleanup_waste','dna_cleanup','end_repair_reagent','blunt','end_repair','a_tail_reagent','tailed','a_tail','adapter','ligation_condition','ligated','ligate','gel_select_condition','gel_select_dna','gel_select_waste','gel_select','forward','reverse','pfu','library','pcr'), g.authority)
    g.add_source_item(SourceItem('unreviewed_original', 'Original supplementary protocols and historical GAII kit not fully reviewed', 'Associated article Methods mentions original one-page cell/tissue supplements and Illumina #1000521/#1000537/#1000538'))
    return g
