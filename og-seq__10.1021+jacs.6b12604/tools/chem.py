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

PAGE_FILE='og-seq.html'
TITLE='OG-seq — oxidative biotin tagging and complementary-strand readout'
NOTES='01_og-seq.html'
SOURCE='<a href="https://doi.org/10.1021/jacs.6b12604">Ding et al. (2017)</a>; associated accepted manuscriptPMC5440228 main24–36. Relevant SupportingInformation remains uninspected.'
CAVEAT='Selected WT-MEF enrichment route. OG-containing and complementary DNA are separate molecules. Tagging creates an oxidized-OG/amine-biotin conjugate, not unchanged OG; exact constitutional chemistry and genomic reaction/workup conditions require SI. No quantitative tagging/capture or perfect purity is predicted. Source sequencing uses complementary DNA released from captured molecules. The historical ssDNA kit, full library and run/index geometry remain unassembled.'
READOUT='Source HiSeq125-cycle paired-end v4, region-scale enrichment near150bp. A sequence read from the complementary strand is not a retained OG chemical base or an exact single-base lesion call. Untreated input control and Ogg1-knockout biological samples are separate source branches.'
PREVIEWS=(('Selected OG enrichment and complementary-strand route','library'),)
@dataclass(frozen=True)
class OGState:
    stage:str
    def __post_init__(self):
        if self.stage not in {'input','tagged','captured','released','bead_retained','library'}:
            raise ValueError('unknown OG-seq molecular stage')
    def rows(self):
        sc=Scene()
        if self.stage in {'input','tagged','captured','bead_retained'}:
            name='original genomic DNA with8-oxoG class; positions unknown' if self.stage=='input' else 'original genomic strand with oxidized-OG/amine-BTN conjugate; exact structure unassembled'
            sc.strand('lesion',[role(name)],label='original lesion-associated strand')
        if self.stage in {'input','tagged','captured','released'}:
            sc.strand('complement',[role('complementary genomic DNA; no covalent BTN tag assigned')],label='complementary strand; sequence unassigned')
        if self.stage=='released':
            sc.strand('carryover',[role('possible other nonbiotin DNA in released fraction; partition unquantified')],label='possible physical carryover, not assumed absent')
        if self.stage=='library':
            sc.strand('library',[role('unassembled source ssDNA-kit end'),role('complement-derived genomic copy'),role('opposite source kit end')],label='source library roles')
            sc.strand('carryover',[role('nonselected/background-derived kit products; contribution unassigned')],label='possible background class')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Associated PMC5440228 main chemical/genomic route24–36; SI/kit boundaries unreviewed')
    dna=g.input('dna','extracted WT-MEF genomicDNA; Ogg1 knockout is a separate biological source',K.MOLECULAR,state=OGState('input'),primary=True)
    shear=g.input('sonication','source sonication near150bp; exact genomic parameters await SI',K.CONDITION)
    fragmented=g.product('fragmented','fragmented DNA before labeling',K.MOLECULAR,state=dna.state)
    g.react('shear','Fragment genomicDNA before chemical tagging',S.FRAGMENTATION,((dna,R.TEMPLATE),(shear,R.CONDITION)),(fragmented,))
    reaction=g.input('reaction','K2IrBr6 and amine-PEG-biotin(BTN); source test5mM/20mM,pH8,75C30min; full genomic recipe unreviewed',K.CONDITION)
    tagged=g.product('tagged','tagged genomic strands and their complements; reagents/workup unassigned',K.MOLECULAR,state=OGState('tagged'))
    g.react('tag','Oxidize OG and trap with amine-BTN',S.CHEMICAL_CONVERSION,((fragmented,R.TEMPLATE),(reaction,R.CONDITION)),(tagged,))
    beads=g.input('beads','source streptavidin-coated magnetic beads; original preparation/washes await SI',K.SELECTION)
    captured=g.product('captured','captured labeled DNA with associated complementary strands',K.MOLECULAR,state=OGState('captured'))
    other=g.product('capture_other','nonretained capture material; exact fraction fate unassembled',K.MOLECULAR,state=molecule('other capture fraction',(role('other DNA/reagents; no perfect partition asserted'),)))
    g.react('capture','Affinity-isolate chemically biotinylated DNA',S.AFFINITY_CAPTURE,((tagged,R.TEMPLATE),(beads,R.SELECTION)),(captured,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Source capture positive; full wash/disposal conditions remain unreviewed.'))
    naoh=g.input('naoh','150mMNaOH20C10min',K.CONDITION)
    released=g.product('released','released complementary DNA for source library preparation',K.MOLECULAR,state=OGState('released'))
    retained=g.product('bead_retained','biotin-adduct strand remains bead-associated',K.MOLECULAR,state=OGState('bead_retained'))
    g.react('release','Release complementary DNA without copying/removing the covalent tagged strand into that fraction',S.STRAND_SEPARATION,((captured,R.TEMPLATE),(naoh,R.CONDITION)),(released,retained))
    g.mark_terminal(Terminal(retained.id,Disposition.UNRESOLVED,'Labeled adduct stable to release conditions; this is not the complementary-strand sequencing substrate. Final disposal awaits SI.'))
    kit=g.input('kit','source single-strandedDNA adaptor kit; identity and exact protocol unreviewed',K.CONDITION)
    library=g.product('library','source complement-derived HiSeq library roles',K.LIBRARY,state=OGState('library'))
    g.react('library_boundary','Historical ssDNA-kit library boundary',S.MAGIC,((released,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Source names ssDNA adaptor preparation and HiSeq endpoint; kit molecular steps, adapters, PCR and read/primer/index assembly are not reconstructed from similarity.')
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source HiSeq PE125 endpoint; complete molecular library and run geometry remain unresolved.'))
    g.declare_reviewed_scope(('dna','sonication','shear','reaction','tag','beads','capture','naoh','release','kit','library_boundary'),g.authority)
    g.add_source_item(SourceItem('supporting_information','Defining published SupportingInformation chemical/genomic recipe and kit identity','NIHMS860150-supplement-Supporting_Information.docx linked from PMC5440228; canonical request gated.'))
    return g
