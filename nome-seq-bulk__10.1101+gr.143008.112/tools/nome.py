"""Source-owned bulk NOMe-seq: explicit methylation calls, physical U, readout T.

No accessibility, enzyme specificity, efficiency or methylation population is
predicted. The reference helper is conditional on caller-supplied methylations.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, workflow_from_protocol_graph, workflow_panel
from page import head, info
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, SequencingPlan, SourceItem)

PAGE_FILE = 'nome-seq.html'
SOURCE = ('Kelly et al.2012 DOI10.1101/gr.143008.112; published PMC3514679 '
          'Results JATS63/75 and Methods119–135; selected IMR90 branch')
READ1_CYCLES = 76


@dataclass(frozen=True)
class MethylatedReference:
    """A selected strand with explicit m5C calls, not a methylation predictor.

    Origins are recorded independently: GCG can have either or both origins.
    Source-described CCG off-target events may also be supplied explicitly.
    """
    dna: str
    endogenous: frozenset[int] = frozenset()
    added: frozenset[int] = frozenset()

    def __post_init__(self):
        if not self.dna or set(self.dna) - set('ACGT'):
            raise ValueError('A reference must be explicit DNA')
        if not isinstance(self.endogenous, frozenset) or not isinstance(self.added, frozenset):
            raise ValueError('Methylation positions must be immutable explicit sets')
        for position in self.endogenous | self.added:
            if type(position) is not int or not 0 <= position < len(self.dna) or self.dna[position] != 'C':
                raise ValueError('Every selected m5C call must identify a reference cytosine')
        for position in self.added:
            if not (position > 0 and self.dna[position - 1] == 'G') and self.context(position) != 'CCG off-target caveat':
                raise ValueError('Only selected source GpC or CCG off-target events are represented')

    def context(self, position):
        if not 0 < position < len(self.dna) - 1 or self.dna[position] != 'C':
            return 'unassigned boundary/context'
        before, after = self.dna[position - 1], self.dna[position + 1]
        if before == 'G':
            return 'GCG ambiguous origin' if after == 'G' else 'GCH'
        if after == 'G':
            return 'CCG off-target caveat' if before == 'C' else 'WCG'
        return 'other context'

    def bisulfite_reference(self):
        protected = self.endogenous | self.added
        return ''.join('U' if base == 'C' and i not in protected else base
                       for i, base in enumerate(self.dna))

    def copied_reference(self):
        """New DNA's reference orientation; not a polymerase-strand orientation claim."""
        return self.bisulfite_reference().replace('U', 'T')


def reference_scene(stage):
    if stage not in {'marked', 'converted', 'copied'}:
        raise ValueError('Unknown conditional base-class stage')
    sc = Scene()
    for ident, label, value in [
        ('unmarked', 'selected unmodified C class', 'U' if stage == 'converted' else 'T' if stage == 'copied' else 'C'),
        ('marked', 'selected m5C class; methylation origin not determined from C alone', 'C')]:
        sc.strand(ident, [Segment('genomic flank', '[genomic flank]', placeholder=True),
                         Segment('physical base', value),
                         Segment('genomic flank', '[genomic flank]', placeholder=True)], label=label)
        if ident == 'marked':
            sc.note(ident, 'Physical base is ' + ('newly copied C' if stage == 'copied' else '5mC') + '; selected idealized relation, not conversion yield.')
    sc.strand('carryover', [Segment('other material', '[reaction mixture]', placeholder=True)],
              label='other molecules remain until source workup; composition unquantified')
    return tuple(sc.rows())


def unknown_ligation():
    sc = Scene()
    sc.strand('adapter', [Segment('historical adapter end', '[unreviewed methylated adapter end]', placeholder=True)])
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,),
                                  'Source specifies methylated Illumina adapters, but historical full polymers/end geometry and linked originals are unreviewed.')


def protocol_graph():
    g = ProtocolGraph('Bulk NOMe-seq — IMR90 source route', authority=SOURCE,
                      representation='published_instance')
    nuclei = g.input('nuclei', 'Isolated IMR90 nuclei; directly cited Miranda2010 isolation protocol unreviewed', K.BIOLOGICAL, primary=True)
    enzyme1 = g.input('enzyme_first', '200U M.CviPI,37°C7.5min; first source dose', K.CATALYST)
    sam1 = g.input('sam_first', 'SAM methyl donor, first source dose; concentration unreviewed', K.MOLECULAR)
    first = g.product('first', 'First-dose chromatin mixture; accessibility/methylation fractions unassigned', K.MOLECULAR)
    g.react('mark_first', 'First methyltransferase incubation', S.CHEMICAL_CONVERSION,
            ((nuclei, R.TEMPLATE), (enzyme1, R.CATALYST), (sam1, R.STRUCTURAL)), (first,))
    enzyme2 = g.input('enzyme_boost', 'Additional100U M.CviPI,37°C7.5min', K.CATALYST)
    sam2 = g.input('sam_boost', 'Fresh SAM methyl donor at boost; concentration unreviewed', K.MOLECULAR)
    marked = g.product('marked', 'Boosted chromatin mixture; both newly marked and unmarked classes persist', K.MOLECULAR,
                       state=reference_scene('marked'))
    g.react('mark_boost', 'Source fresh enzyme/SAM boost', S.CHEMICAL_CONVERSION,
            ((first, R.TEMPLATE), (enzyme2, R.CATALYST), (sam2, R.STRUCTURAL)), (marked,))
    quench = g.input('stop_reagent', 'Source reaction stop; exact reagent/procedure unreviewed', K.CONDITION)
    stopped = g.product('stopped', 'Stopped source nuclear methylation reaction', K.MOLECULAR, state=marked.state)
    g.react('stop', 'Stop reaction before extracting DNA', S.QUENCHING,
            ((marked, R.TEMPLATE), (quench, R.CONDITION)), (stopped,))
    extraction = g.input('extraction', 'Source DNA extraction; exact historical procedure unreviewed', K.SELECTION)
    dna = g.product('dna', 'Extracted DNA,5µg used for library construction', K.MOLECULAR, state=marked.state)
    waste = g.product('extraction_waste', 'Physically rejected extraction fractions; exact composition unassigned', K.WASTE)
    g.react('extract', 'Extract source marked genomic DNA', S.EXTRACTION,
            ((stopped, R.TEMPLATE), (extraction, R.SELECTION)), (dna, waste))
    g.discard(waste, 'Physical extraction fractions; no specific quantitative depletion claim')
    fragmenter = g.input('fragmentation', 'Source fragmentation to approximately200bp; exact procedure unreviewed', K.CONDITION)
    fragments = g.product('fragments', 'Genomic fragment mixture; no selection inferred from target size', K.MOLECULAR, state=dna.state)
    g.react('fragment', 'Fragment extracted DNA before library end processing', S.FRAGMENTATION,
            ((dna, R.TEMPLATE), (fragmenter, R.CONDITION)), (fragments,))
    repair = g.input('repair_enzyme', 'Epicentre end-repair reagent; precise historical formulation unreviewed', K.CATALYST)
    repaired = g.product('repaired', 'Source end-repaired DNA; exact terminal geometry unassembled', K.MOLECULAR, state=dna.state)
    g.react('repair', 'End-repair source fragments', S.END_REPAIR,
            ((fragments, R.TEMPLATE), (repair, R.CATALYST)), (repaired,))
    adapters = g.input('adapters', 'Historical methylated Illumina adapters; polymers/modifications unassembled', K.ADAPTER)
    ligase = g.input('ligation_reagents', 'Source ligation formulation; historical linked-original details unreviewed', K.CONDITION)
    ligated = g.product('ligated', 'Methylated adapter-bearing genomic DNA; possible other molecular classes unquantified', K.MOLECULAR, state=dna.state)
    g.react('ligate', 'Ligate methylated adapters BEFORE bisulfite conversion', S.LIGATION,
            ((repaired, R.TEMPLATE), (adapters, R.STRUCTURAL), (ligase, R.CONDITION)), (ligated,), geometry=unknown_ligation())
    bisulfite = g.input('bisulfite', 'Zymo EZ DNA Methylation source procedure', K.CONDITION)
    converted = g.product('converted', 'Bisulfite-treated DNA; selected unmodified-C→U reference class', K.MOLECULAR, state=reference_scene('converted'))
    g.react('convert', 'Bisulfite conversion before PCR; physical U distinct from reported T', S.CHEMICAL_CONVERSION,
            ((ligated, R.TEMPLATE), (bisulfite, R.CONDITION)), (converted,))
    workup = g.input('kit_workup', 'Historical Zymo kit workup and pre-PCR recovery not read', K.CONDITION)
    ready = g.product('ready', 'Source bisulfite material supplied to PCR; exact carryover/ends unassigned', K.MOLECULAR, state=converted.state)
    g.react('workup', 'Unassembled historical kit workup boundary', S.MAGIC,
            ((converted, R.TEMPLATE), (workup, R.CONDITION)), (ready,),
            reason='Source names the bisulfite kit; no undocumented wash, debris disposal, end repair or final purity is invented.')
    f = g.input('forward', 'Historical source PCR forward primer; polymer unreviewed', K.OLIGO)
    r = g.input('reverse', 'Historical source PCR reverse primer; polymer unreviewed', K.OLIGO)
    pcr = g.input('pcr_condition', 'Source PCR formulation/program not inspected in directly linked originals', K.CONDITION)
    library = g.product('library', 'PCR-amplified NOMe library; copies and original template/carryover not quantitatively resolved', K.LIBRARY,
                        state=reference_scene('copied'))
    g.react('pcr', 'PCR copies U-derived sequence as T; methylated-origin bases copy as C', S.PCR,
            ((ready, R.TEMPLATE), (f, R.STRUCTURAL), (r, R.STRUCTURAL), (pcr, R.CONDITION)), (library,),
            geometry=PCRGeometry.unavailable(Scene(), reason='Historical full template/primer polymers, physical strand assembly and intervening workup unassembled.'))
    g.sequence(library, SequencingPlan(acquired_roles=('Read 1',), read_lengths={'Read 1':READ1_CYCLES},
                                     limitation='IMR90 source HiSeq2000 single-end76; complete library/run-primer/index geometry unassembled. Separate GBM paired-end50/100 runs are not merged into this plan.'))
    g.declare_reviewed_scope(('nuclei','enzyme_first','sam_first','mark_first','enzyme_boost','sam_boost','mark_boost','stop_reagent','stop',
                             'extraction','extract','fragmentation','fragment','repair_enzyme','repair','adapters','ligation_reagents',
                             'ligate','bisulfite','convert','kit_workup','workup','forward','reverse','pcr_condition','pcr'), SOURCE)
    g.add_source_item(SourceItem('linked_source_boundaries', 'Nuclear isolation, complete historical adapter/end/library/run geometry and supplemental originals',
                                'Directly cited Miranda2010/Lister2009/Kelly2010/Berman2012 and publisher supplement index remain unreviewed or blocked.'))
    return g


def render_page():
    return ('<!doctype html><html lang="en">' + head('NOMe-seq — source-limited bulk model') +
            '<main class="wrap"><h1>NOMe-seq — bulk IMR90 branch</h1>' +
            '<p><a href="01_nome-seq.html">Research notes</a></p>' +
            info('Private incomplete model. Conditional methylation/base classes, not a prediction of nucleosome occupancy or conversion efficiency. '
                 'GCG origin is ambiguous and source-reported CCG off-target marking is not excluded. '
                 'Historical end/library geometry and linked originals remain unassembled.') +
            workflow_panel(workflow_from_protocol_graph(protocol_graph(), 'library')) +
            '<p>GCH accessibility, HCG/CCG caveats, SNP checks, duplicate removal and GCG exclusion are computational interpretation; they do not remove physical bases. '
            'Source IMR90 uses single-end76. Glioblastoma libraries have separate paired-end50/100 schedules.</p></main></html>')


def selftest():
    from checks import Check
    c = Check()
    # Synthetic calculation fixtures, not study molecules or source oligos.
    reference = MethylatedReference('AGCAGCGACGACCG', frozenset({8}), frozenset({2,5,12}))
    c('context distinctions retain GCG ambiguity and CCG off-target caveat',
      tuple(reference.context(i) for i in (2,5,8,12)),
      ('GCH','GCG ambiguous origin','WCG','CCG off-target caveat'))
    c('physical bisulfite product has U before new DNA reports T',
      (reference.bisulfite_reference()[11], reference.copied_reference()[11]), ('U','T'))
    c('new DNA sequence preserves length and records supplied protected positions',
      (len(reference.copied_reference()), tuple(reference.copied_reference()[i] for i in (2,5,8,12))),
      (len(reference.dna), ('C','C','C','C')))
    edge = MethylatedReference('GC', added=frozenset({1}))
    c('a supplied terminal GpC mark does not invent its unavailable third context base',
      (edge.copied_reference(), edge.context(1)), ('GC','unassigned boundary/context'))
    for label, args in [('non-C methylation', ('ACG',frozenset({0}),frozenset())),
                        ('unsupported added site', ('ACA',frozenset(),frozenset({1}))),
                        ('out-of-range site', ('GCA',frozenset({9}),frozenset()))]:
        try:
            MethylatedReference(*args)
        except ValueError:
            rejected=True
        else:
            rejected=False
        c('selected-call boundary rejects '+label,rejected,True)
    c('unreviewed originals prevent finished publication',protocol_graph().publishable,False)
    c.report()


if __name__ == '__main__':
    selftest()
