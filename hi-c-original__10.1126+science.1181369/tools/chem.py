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
    c('Actual unread SI obligation stays uncovered','pending_assembly' in g.source_coverage,False)
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c('Every physical DNA output has a declared downstream operation or terminal', not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c.report()

PAGE_FILE='hi-c-original.html'
TITLE='Hi-C (original2009)'
NOTES='01_hi-c-original.html'
SOURCE='<a href="https://doi.org/10.1126/science.1181369">Defining article</a>; Associated author manuscript PMC2858594, main112–118 and Fig1caption320'
CAVEAT='Main-only selected HindIII branch. Biotin nucleotide identity, enzyme/end geometry, cleanup recipe and final adapters require the linked gated Supplementary text. The NcoI replicate is separate. No in-situ or DNase Hi-C recipe is inherited.'
READOUT='Illumina Genome Analyzer paired-end contact-pair sequencing; full adapter/run/index geometry unassembled.'
PREVIEWS=(("Source-qualified partial molecular route","output"),)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    current=g.input('dna','source genomic DNA/chromatin; genomic bases unassigned',K.MOLECULAR,state=molecule('source DNA',(role('genomic DNA; coordinates and sequence unassigned'),)),primary=True)
    reagent=g.input('fix_reagent','Formaldehyde crosslink spatially adjacent chromatin: actual molecular reagent; source dose unassigned',K.MOLECULAR)
    product=g.product('fix_product','Formaldehyde crosslink spatially adjacent chromatin',K.MOLECULAR,state=molecule('source product classes',(role('crosslinked DNA and associated protein complexes; no fixed genomic coordinates'),)))
    fix_conditions=g.input('fix_conditions','Source formaldehyde fixation; dose/time/temperature remain unassigned pending actual SI, no invented conditions',K.CONDITION)
    g.react('fix','Formaldehyde crosslink spatially adjacent chromatin',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(reagent,R.STRUCTURAL),(fix_conditions,R.CONDITION)),(product,))
    current=product
    reagent=g.input('digest_reagent','HindIII digestion; source five-prime-overhang class',K.CATALYST)
    product=g.product('digest_product','HindIII digestion; source five-prime-overhang class',K.MOLECULAR,state=molecule('source product classes',(role('crosslinked restriction fragments with source5′overhangs; small fragments retained'),)))
    g.react('digest','HindIII digestion; source five-prime-overhang class',S.DIGESTION,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,))
    current=product
    reagent=g.input('fill_reagent','Actual fill-in nucleotide-containing supply including biotin residue; nucleotide identity and enzyme inventory pending genuine SI',K.MOLECULAR)
    product=g.product('fill_product','Fill overhangs including a biotinylated residue; nucleotide identity pending SI',K.MOLECULAR,state=molecule('source product classes',(role('blunt genomic ends with incorporated biotinylated nucleotide; exact position and molecular assembly unreviewed'),)))
    g.react('fill','Fill overhangs including a biotinylated residue; nucleotide identity pending SI',S.MAGIC,((current,R.TEMPLATE),(reagent,R.STRUCTURAL)),(product,),reason='Main specifies fill-in; biotin nucleotide identity and exact molecular geometry await SI.')
    current=product
    reagent=g.input('join_reagent','Source dilute blunt-end ligation',K.CATALYST)
    product=g.product('join_product','Source dilute blunt-end ligation',K.MOLECULAR,state=molecule('source product classes',(role('ligation-product population including biotin-marked contact junctions and unligated/background classes'),)))
    g.react('join','Source dilute blunt-end ligation',S.PROXIMITY_LIGATION,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,))
    current=product
    current=cleanup(g,'purify_product','Source DNA purification, recipe pending SI',current,molecule('recovered DNA',(role('purified ligation-product population; exact partition unquantified'),)))
    reagent=g.input('shear_reagent','Shear purified DNA',K.CONDITION)
    product=g.product('shear_product','Shear purified DNA',K.MOLECULAR,state=molecule('source product classes',(role('sheared junction-containing and nonjunction DNA fragments'),)))
    g.react('shear','Shear purified DNA',S.FRAGMENTATION,((current,R.TEMPLATE),(reagent,R.CONDITION)),(product,))
    current=product
    reagent=g.input('capture_reagent','Source streptavidin selection',K.SELECTION)
    product=g.product('output','Source streptavidin selection',K.MOLECULAR,state=molecule('source product classes',(role('streptavidin-selected biotin-containing DNA; no perfect contact-product specificity assumed'),)))
    other=g.product('capture_other','physically separated source flowthrough/washes; partition unquantified',K.WASTE)
    g.react('capture','Source streptavidin selection',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(reagent,R.SELECTION)),(product,other))
    g.discard(other,'Actual source isolation/washes; no perfect specificity asserted.')
    current=product
    g.mark_terminal(Terminal(current.id,Disposition.UNRESOLVED,'Supplementary text and Supplementary Figures linked from PMC2858594, canonical requests returned PMC preparing-download gates.'))
    # Literal main-only operational inventory; this does not cover unavailable SI.
    g.declare_reviewed_scope(('dna','fix_reagent','fix_conditions','fix','digest_reagent','digest','fill_reagent','fill','join_reagent','join','purify_product_condition','purify_product_step','shear_reagent','shear','capture_reagent','capture'),g.authority)
    g.add_source_item(SourceItem('pending_assembly','Unreviewed source/assembly boundary','Supplementary text and Supplementary Figures linked from PMC2858594, canonical requests returned PMC preparing-download gates.'))
    return g
