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
    c('Every DNA molecular output continues or has explicit terminal fate',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    order=list(g.steps)
    c('Pre-existing end blocking precedes nuclease and new-end biotin tailing',order.index('block_original_ends')<order.index('mungbean')<order.index('biotin_tail'))
    c('Nuclease action retains physical DNA products',all(g.materials[mid].kind is not K.WASTE for mid in g.steps['mungbean'].outputs))
    c('Fragmentation retains both labelled and unlabelled source classes',all(g.materials[mid].kind is not K.WASTE for mid in g.steps['sonicate'].outputs))
    c('Capture companion fractions have unresolved source fate rather than invented disposal',g.terminals['nonselected'].disposition,Disposition.UNRESOLVED)
    c.report()

def population(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other template strands, RNA/protein or reaction fragments; actual molecular recovery unquantified'),)))

def op(g,ident,label,kind,current,new_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,*,additional=()):
    reagent=g.input(ident+'_reagent',reagent_label,reagent_kind)
    product=g.product(ident+'_product',label,K.MOLECULAR,state=new_state)
    kw={'reason':'Source operation established; exact molecular reaction geometry and population fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((current,R.TEMPLATE),(reagent,reagent_role))+tuple(additional),(product,),**kw)
    return product

PAGE_FILE='ssdna-seq.html'
TITLE='ssDNA-Seq (permanganate)'
NOTES='01_ssdna-seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.cell.2013.04.033">Defining article</a>; Associated accepted manuscript PMC3684982 Results104–110/Fig3 and Methods146–148; detailed originalSI gated'
CAVEAT='Selected original2013 permanganate/mung-bean/TdT route. Oxidation, blocked original3′ends and new nuclease-associated/TdT-accepted ends are separate source classes; exact lesion constitution,5′geometry, quantities and workups require originalSI. Later2017S1/detailed library recipe is not silently substituted. Molecular digestion products remain until actual source separation. INFERRED — TdT-accepted3′OH is a selected end interpretation from source tailing; exact nuclease end geometry remains unreviewed.'
READOUT='Source Illumina deep sequencing of biotin-selected fragments near200bp; exact adapter/library and acquired read/cycle/index geometry remain unassembled.'
PREVIEWS=(("Source-qualified molecular operation route","library"),)

@dataclass(frozen=True)
class EndHistory:
    stage:str
    def __post_init__(self):
        if self.stage not in {'native','oxidized','blocked','cleaved','tagged','sheared','captured','library'}:raise ValueError('unsupported original permanganate ssDNA-seq stage')
    def rows(self):
        sc=Scene()
        if self.stage=='library':
            sc.strand('library',[role('unassembled source Illumina library end'),role('genomic end-associated DNA fragment'),role('opposite source end; biotin-tail processing unreviewed')],label='library roles, no complete sequence claim')
        else:
            label='source genomic DNA; no oxidation assigned' if self.stage=='native' else 'source oxidized pyrimidine class, preferentially exposedT; exact adduct constitution unassembled'
            sc.strand('genomic',[role(label)],label='original genomic DNA class')
            if self.stage in {'blocked','cleaved','tagged','sheared','captured'}:
                sc.strand('old_end',[role('selected pre-existing3′end chain-terminated by source cordycepin/TdT; exact appended nucleotide geometry unassembled')],label='blocked original-end class, not5′phosphorylation')
            if self.stage in {'cleaved','tagged','sheared','captured'}:
                end='new nuclease-associated DNA end; TdT-accepted3′OH class is a reaction-derived interpretation, exact5′geometry unreviewed'
                if self.stage in {'tagged','sheared','captured'}:end='new DNA end with source3′TdT biotin-nucleotide tail; exact sequence/length/nucleotide identity unassembled'
                sc.strand('new_end',[role(end,inferred=self.stage=='cleaved')],label='new end class;5′end not assigned')
                sc.strand('digestion_products',[role('other genomic fragments and possible cleavage products; no disposal at nuclease action')],label='other physical DNA')
        sc.strand('carryover',[role('possible unblocked/unlabelled/backgroundDNA, RNA/protein and reaction remnants; fractions unquantified')],label='possible other material, no perfect specificity')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    current=g.input('cells','source Raji or mouseBcells; separate biological experiments, not an implicit physical pool',K.MOLECULAR,state=EndHistory('native'),primary=True)
    exposure=g.input('permanganate_environment','Source living-cell oxidation environment; actual concentration, time and temperature unreviewed in original SI',K.CONDITION)
    current=op(g,'permanganate','Source living-cell KMnO4 oxidation, preferential exposed pyrimidines/T; no arbitrary genomic coordinate assigned',S.CHEMICAL_CONVERSION,current,EndHistory('oxidized'),'Actual source KMnO4 molecular oxidant; dose/workup and exact originalSI recipe uninspected',K.MOLECULAR,R.STRUCTURAL,additional=((exposure,R.CONDITION),))
    current=op(g,'recovery_boundary','Original genomicDNA isolation/quenches/workups before end processing; detailedSI uninspected',S.MAGIC,current,current.state,'source original preparation boundary, not2017S1recipe; physical partition and RNase/protease details unreviewed')
    tdt=g.input('blocking_tdt','source TdT for original-end chain termination',K.CATALYST)
    cord=g.input('cordycepin','source chain-terminating cordycepin5′triphosphate',K.MOLECULAR,state=molecule('source blocking nucleotide',(role('cordycepin5′triphosphate; exact incorporation assembly unreviewed'),)))
    blocked=g.product('blocked','selected blocked original3′DNAends plus background/unblocked material',K.MOLECULAR,state=EndHistory('blocked'))
    g.react('block_original_ends','Source cordycepin/TdT blocks pre-existing3′ends BEFORE mungbean cleavage',S.MAGIC,((current,R.TEMPLATE),(tdt,R.CATALYST),(cord,R.STRUCTURAL)),(blocked,),reason='Source chain-termination event established; exact end chemistry/efficiency and intervening cleanup from gated originalSI unassembled.')
    current=op(g,'mungbean','INFERRED — selected TdT-accepted3OH class after source mungbean cleavage; exact complete termini unreviewed',S.DIGESTION,blocked,EndHistory('cleaved'),'source mungbean nuclease; exact5′end polarity/bases/conditions unreviewed',K.CATALYST,R.CATALYST)
    tdt=g.input('tagging_tdt','source TdT for newly exposed3′end tailing',K.CATALYST)
    nucleotides=g.input('biotin_nucleotides','source biotinylated nucleotide tailing pool; exact original nucleotides/ratios unreviewed',K.MOLECULAR,state=molecule('source tailing nucleotides',(role('biotinylated nucleotides; no index or UMI meaning assigned'),)))
    tagged=g.product('tagged','source3′biotin-tailed end-associated DNA plus unlabelled/other fragments',K.MOLECULAR,state=EndHistory('tagged'))
    g.react('biotin_tail','Source TdT3′tailing of nuclease-associated ends, not covalent tag transfer to every genomic molecule',S.TAILING,((current,R.TEMPLATE),(tdt,R.CATALYST),(nucleotides,R.STRUCTURAL)),(tagged,))
    current=op(g,'sonicate','Source sonication near200bp; labelled and unlabelled DNA remain',S.FRAGMENTATION,tagged,EndHistory('sheared'),'source sonication; exact conditions requireSI')
    beads=g.input('streptavidin','source streptavidin selection of biotinylated DNA',K.SELECTION)
    selected=g.product('captured','source biotin-selected DNA class; capture/recovery and background unquantified',K.MOLECULAR,state=EndHistory('captured'))
    other=g.product('nonselected','source nonselected companion fractions; exact physical fate/workup unassembled',K.MOLECULAR,state=population('unlabelled/other genomic fragments and preparation material'))
    g.react('capture','Source physical streptavidin selection after sonication',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(beads,R.SELECTION)),(selected,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Do not assign uninspected original wash/disposal recipe.'))
    prep=g.input('illumina_prep','Actual source-delegated Illumina library molecular supplies; stock packaging, tail removal/end repair/adapters/PCR unreviewed',K.MOLECULAR)
    library=g.product('library','source Illumina end-associated genomic library roles',K.LIBRARY,state=EndHistory('library'))
    g.react('library_boundary','Original genomic Illumina library boundary; no2017S1/tail-removal/adapter recipe inherited',S.MAGIC,((selected,R.TEMPLATE),(prep,R.STRUCTURAL)),(library,),reason='Source main/Fig3 establishes Illumina endpoint; gated2013SI determines actual molecular construction and biotin-tail fate.')
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source deep sequencing; full strand/end/run geometry remains unassembled.'))
    # Literal main-supported inputs/operations; unread workups and library
    # construction are deliberately outside this source inventory.
    g.declare_reviewed_scope(('cells','permanganate_reagent','permanganate_environment','permanganate','blocking_tdt','cordycepin','block_original_ends','mungbean_reagent','mungbean','tagging_tdt','biotin_nucleotides','biotin_tail','sonicate_reagent','sonicate','streptavidin','capture'),'Associated manuscript PMC3684982 Results104–106. Selected operational classes only; exact molecular ends, doses, recovery and Illumina producer remain unreviewed.')
    g.add_source_item(SourceItem('gated_original_si','Original2013ExtendedExperimentalProcedures/figures and tables','https://pmc.ncbi.nlm.nih.gov/articles/instance/3684982/bin/NIHMS472536-supplement-01.pdf ; actual preparing-download gate. Companion02.xls also gated.'))
    return g
