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
    before=G4Observation('ACGGGT','Na+','ACGGGT');after=before.under('Na++PDS','ACNNNN')
    c('Condition-dependent observed calls cannot overwrite original DNA sequence',(after.template_dna,before.template_dna),('ACGGGT','ACGGGT'))
    c('Changed read calls remain separate measured values',after.reported_calls,'ACNNNN')
    c('Both source stabilized-condition graphs remain incomplete',all(graph_for_condition(v).publication_problems() for v in ('K+','Na++PDS')),True)
    c.report()

def population(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other template strands, RNA/protein or reaction fragments; actual molecular recovery unquantified'),)))

def op(g,ident,label,kind,current,new_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    reagent=g.input(ident+'_reagent',reagent_label,reagent_kind)
    product=g.product(ident+'_product',label,K.MOLECULAR,state=new_state)
    kw={'reason':'Source operation established; exact molecular reaction geometry and population fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((current,R.TEMPLATE),(reagent,reagent_role)),(product,),**kw)
    return product

PAGE_FILE='g4-seq.html'
TITLE='G4-seq'
NOTES='01_g4-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nbt.3295">Defining article</a>; Published publisher title/abstract; published SI MOESM27Fig1–7/captions,MOESM28Tables1–3 and experimental-result scope'
CAVEAT='Selected original2015 DNA-sequencing interference assay; not RNA reverse transcription or2019G4-seq2. Full main Methods/library construction remain blocked. Sequencing scripts31–33 were fetched but unreviewed under cap; exact repeat-read/strand/index/run geometry is not assigned from standard paired-end conventions. Conditions change DNA conformation and sequencing observations, not the original DNA base identity.'
READOUT='Source controls compare first-readNa+ with second-readNa++PDS orK+-stabilized conditions; exact full physical acquired read geometry remains unassembled.'
PREVIEWS=(("Source-qualified molecular operation route","library"),)

@dataclass(frozen=True)
class G4Observation:
    template_dna:str
    condition:str
    reported_calls:str
    def __post_init__(self):
        if not self.template_dna or set(self.template_dna)-set('ACGT') or self.condition not in {'Na+','K+','Na++PDS'} or set(self.reported_calls)-set('ACGTN'):raise ValueError('explicit DNA, source condition and caller-supplied observed calls required')
    def under(self,condition,calls):
        return G4Observation(self.template_dna,condition,calls)

@dataclass(frozen=True)
class SequencingPopulation:
    stage:str
    def rows(self):
        if self.stage not in {'input','library','cluster','baseline','stabilized'}:raise ValueError('unsupported G4-seq source stage')
        sc=Scene()
        sc.strand('template',[role('original genomic DNA-derived template sequence; never converted to erroneous base calls')],label='DNA template, unassigned sequence/ends')
        if self.stage in {'baseline','stabilized'}:sc.strand('copies',[role('nascent sequencing-extension strands; exact orientation/length and stall classes unassembled')],label='sequencing products, not altered original template')
        sc.strand('other',[role('other library strands, reagents and possible unextended/stalled DNA; exact population unassigned')],label='other physical molecules')
        return sc.rows()

def graph_for_condition(condition='Na++PDS'):
    if condition not in {'K+','Na++PDS'}:raise ValueError('two selected original stabilized-condition branches only')
    g=ProtocolGraph(TITLE+' — '+condition,authority=SOURCE,representation='published_instance')
    current=g.input('genomic_dna','source human genomic DNA; full sample/preparation details require main Methods',K.MOLECULAR,state=SequencingPopulation('input'),primary=True)
    current=op(g,'library_boundary','Source genomic sequencing-library preparation; main Methods uninspected',S.MAGIC,current,SequencingPopulation('library'),'original adapter/end/PCR/selection operations not reconstructed from matching Illumina controls')
    current=op(g,'flowcell_boundary','Source independent sequencing clusters; exact bridge/strand/run geometry unassembled',S.MAGIC,current,SequencingPopulation('cluster'),'source Illumina array-sequencing preparation and original scripts; not a standard paired-end default')
    current=op(g,'first_read','Source first-readNa+ sequencing; observed reference-like calls do not mutate template',S.MAGIC,current,SequencingPopulation('baseline'),'source Na+ sequencing buffer/polymerase/primer system, exact chemistry unreviewed')
    current=op(g,'stabilize','Source change to '+condition+' G4-stabilizing sequencing conditions; no guaranteed fold or DNA base conversion',S.GENERIC,current,SequencingPopulation('stabilized'),'source '+condition+' buffer/ligand; exact exchange/denaturation/strand cycle program unreviewed')
    final=g.product('library','Source second-read sequencing population; original genomic DNA template retained',K.LIBRARY,state=SequencingPopulation('stabilized'))
    system=g.input('second_read_system','source second-read sequencing reaction; exact reset/run geometry requires original scripts',K.CONDITION)
    g.react('second_read','Source stabilized second-read DNApolymerase interference; mismatch/quality is measured data, not a template mutation',S.MAGIC,((current,R.TEMPLATE),(system,R.CONDITION)),(final,),reason='Source SI Fig3–5 establishes compared reads and suspected polymerase stalling; exact reset/copy/primer/run geometry unreviewed.')
    g.mark_terminal(Terminal(final.id,Disposition.UNRESOLVED,'Source G4-seq genomic array readout; do not assign a standard paired-end orientation or complete library.'))
    g.add_source_item(SourceItem('blocked_main','Defining full main Methods, complete sequencing-library construction','https://doi.org/10.1038/nbt.3295 ; standard fetch reports no open fulltext.'))
    g.add_source_item(SourceItem('unreviewed_scripts','Original sequencing XML scripts and complete physical read geometry','Published MOESM31/32/33 exact originals fetched but unreviewed within cap; scripts differ by instrument/condition.'))
    return g

def protocol_graph():return graph_for_condition()
