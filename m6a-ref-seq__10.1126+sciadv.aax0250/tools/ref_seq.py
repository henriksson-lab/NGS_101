"""Published m6A-REF-seq matched FTO/MazF chemistry before historical NR801 kit.

This owner does not inherit MAZTER-seq's adapters, RT or library chemistry.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, oligo, workflow_from_protocol_graph, workflow_panel
from page import head, info
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, Terminal, Disposition, SourceItem)

PAGE_FILE = 'm6a-ref-seq.html'
SOURCE = ('Zhang et al.2019, DOI10.1126/sciadv.aax0250; published Methods '
          'RNA materials and m6A-REF-seq library construction, PMC6609220 JATS139–169')
CONTROL_RNA = 'UUGGUUUUUUUUGGACAUGUAUAUAGU'
PCR_CYCLES = 15


@dataclass(frozen=True)
class RNAHistory:
    """History of a population, not conversion, cleavage or methylation efficiency."""
    fto_exposed: bool = False
    mazf_exposed: bool = False
    end_repaired: bool = False

    def __post_init__(self):
        if self.end_repaired and not self.mazf_exposed:
            raise ValueError('End-repaired endpoint belongs after MazF in this represented route')

    def demethylation_treatment(self):
        if self.mazf_exposed or self.fto_exposed:
            raise ValueError('The represented FTO pretreatment precedes MazF, once')
        return replace(self, fto_exposed=True)

    def mazf_treatment(self):
        if self.mazf_exposed:
            raise ValueError('The represented route has one MazF treatment')
        return replace(self, mazf_exposed=True)

    def repair(self):
        if not self.mazf_exposed or self.end_repaired:
            raise ValueError('The represented end repair follows MazF, once')
        return replace(self, end_repaired=True)

    def rows(self):
        sc = Scene()
        if self.mazf_exposed:
            sc.strand('uncut', [Segment('retained uncut RNA', '[uncut RNA]', placeholder=True)],
                      label='uncut/protected molecules remain')
            sc.strand('left', [Segment('upstream RNA fragment', '[upstream RNA]', placeholder=True)],
                      label='upstream product; terminal chemistry unassembled')
            sc.strand('right', [Segment('recognition motif', 'ACA'),
                                Segment('downstream RNA fragment', '[downstream RNA]', placeholder=True)],
                      label='selected downstream product; terminal chemistry unassembled')
            sc.note('right', 'Selected ACA event begins at the first A; not an exhaustive or quantitative digest prediction.')
        else:
            sc.strand('rna', [Segment('RNA population', '[RNA population]', placeholder=True)],
                      label='RNA modification states and fractions unassigned')
        first = 'uncut' if self.mazf_exposed else 'rna'
        if self.fto_exposed:
            sc.note(first, 'FTO-treated population; residual m6A/intermediates and conversion fractions unassigned, not assumed fully erased.')
        if self.end_repaired:
            sc.note(first, 'Source T4PNK treatment recorded; exact phosphate/hydroxyl state remains unassembled.')
        return tuple(sc.rows())


def control_oligo():
    """Verbatim unmodified validation control; not inserted into the NGS branch."""
    return oligo('RNA-ACA validation substrate', [Segment('source RNA-ACA', CONTROL_RNA)])


def protocol_graph():
    g = ProtocolGraph('m6A-REF-seq — matched FTO and buffer-control prefixes', authority=SOURCE,
                      representation='published_instance')
    original = RNAHistory()
    mrna = g.input('mrna', 'Source poly(A) RNA; HEK293T Dynabeads purification or purchased tissue mRNA',
                   K.MOLECULAR, primary=True, state=original.rows(), chemical_state=original)
    preheat = g.input('initial_heat', '85°C for 5 min, then ice for 2 min, before FTO pretreatment', K.CONDITION)
    initial = g.product('initial_denatured', 'Initially heat-treated mRNA; post-cooling structure unassigned', K.MOLECULAR,
                        state=original.rows(), chemical_state=original)
    g.react('initial_denature', 'First source heat treatment, before FTO or matched-buffer preparation', S.THERMAL_TREATMENT,
            ((mrna, R.TEMPLATE), (preheat, R.CONDITION)), (initial,))
    control = g.product('control', 'Matched minus-FTO aliquot', K.MOLECULAR, state=original.rows(), chemical_state=original)
    treatment = g.product('treatment', 'Matched plus-FTO aliquot', K.MOLECULAR, state=original.rows(), chemical_state=original)
    g.react('split', 'Allocate source matched FTO and minus-enzyme controls', S.SPLIT,
            ((initial, R.TEMPLATE),), (control, treatment))
    reviewed = ['mrna', 'initial_heat', 'initial_denature', 'split']
    for name, template, exposed in [('minus', control, False), ('plus', treatment, True)]:
        buffer = g.input(name + '_buffer', 'Source Fe(II), alpha-ketoglutarate, ascorbate, RNase inhibitor and Tris reaction buffer', K.CONDITION)
        inputs = [(template, R.TEMPLATE), (buffer, R.CONDITION)]
        history = original
        if exposed:
            enzyme = g.input(name + '_fto', '2.5 µg FTO demethylase; source 100–200 ng mRNA in 20 µL', K.CATALYST)
            inputs.append((enzyme, R.CATALYST))
            history = history.demethylation_treatment()
        treated = g.product(name + '_treated', 'FTO-treated RNA' if exposed else 'Minus-FTO matched-buffer RNA',
                            K.MOLECULAR, state=history.rows(), chemical_state=history)
        g.react(name + '_pretreat', 'Source 3 h at room temperature; EDTA quench' if exposed else 'Source matched reaction mixture without FTO protein',
                S.CHEMICAL_CONVERSION if exposed else S.THERMAL_TREATMENT, inputs, (treated,))
        heat = g.input(name + '_heat', '85°C for 5 min, ice for 2 min', K.CONDITION)
        heated = g.product(name + '_heated', 'Heat-treated RNA and carried pretreatment mixture', K.MOLECULAR,
                           state=history.rows(), chemical_state=history)
        g.react(name + '_denature', 'Heat both matched samples before MazF', S.THERMAL_TREATMENT,
                ((treated, R.TEMPLATE), (heat, R.CONDITION)), (heated,))
        mazf = g.input(name + '_mazf', '10 U MazF in source buffer, 37°C for 30 min', K.CATALYST)
        history = history.mazf_treatment()
        digested = g.product(name + '_digested', 'RNA fragments plus uncut/protected RNA; no computational filtering yet',
                             K.MOLECULAR, state=history.rows(), chemical_state=history)
        g.react(name + '_digest', 'MazF digest; selected ACA products and uncut molecules represented', S.DIGESTION,
                ((heated, R.TEMPLATE), (mazf, R.CATALYST)), (digested,))

        def clean(identifier, material, state, label):
            purified = g.product(identifier, label, K.MOLECULAR, state=state.rows(), chemical_state=state)
            waste = g.product(identifier + '_waste', 'Physical RNA Clean & Concentrator-5 flowthrough/washes', K.WASTE)
            g.react(identifier + '_cleanup', 'Source Zymo RNA Clean & Concentrator-5 purification', S.CLEANUP,
                    ((material, R.TEMPLATE),), (purified, waste))
            g.discard(waste, 'Physical purification; exact population partition and quantitative recovery not asserted')
            return purified

        purified = clean(name + '_clean', digested, history, 'Purified fragmented/uncut RNA population')
        pnk = g.input(name + '_pnk', 'Vazyme T4PNK, source 1x buffer; 37°C for 30 min', K.CATALYST)
        history = history.repair()
        repaired = g.product(name + '_repaired', 'Source end-repaired RNA; precise termini unassigned', K.MOLECULAR,
                             state=history.rows(), chemical_state=history)
        g.react(name + '_repair', 'Source T4PNK end repair', S.END_REPAIR,
                ((purified, R.TEMPLATE), (pnk, R.CATALYST)), (repaired,))
        ready = clean(name + '_ready', repaired, history, 'Purified source RNA ready for historical VAHTS NR801 kit')
        g.mark_terminal(Terminal(ready.id, Disposition.UNRESOLVED,
                                f'Source uses 10ng RNA in VAHTS NR801 small-RNA kit, PCR{PCR_CYCLES}cycles,8%urea PAGE then HiSeqX10; historical kit polymers/reaction geometry and run configuration unassembled.'))
        reviewed.extend(name + suffix for suffix in ['_buffer', '_pretreat', '_heat', '_denature', '_mazf', '_digest', '_clean_cleanup', '_pnk', '_repair', '_ready_cleanup'])
        if exposed:
            reviewed.append(name + '_fto')
    g.declare_reviewed_scope(reviewed, SOURCE)
    g.add_source_item(SourceItem('historical_nr801', 'Historical VAHTS NR801 kit molecular steps, oligos, PCR and run geometry',
                                'Published library-construction Methods name the kit; its relevant historical original is not yet identified/inspected.'))
    return g


def render_page():
    g = protocol_graph()
    return ('<!doctype html><html lang="en">' + head('m6A-REF-seq — source-limited prefixes') +
            '<main class="wrap"><h1>m6A-REF-seq</h1>' +
            '<p><a href="01_m6a-ref-seq.html">Research notes</a></p>' +
            info('Private incomplete model. FTO and matched minus-enzyme pretreatments remain separate. '
                 'Population conversion/cleavage is not predicted. Historical NR801 kit and run geometry are unassembled; '
                 'MAZTER-seq adapters/RT are not imported.') +
            '<h2>Minus FTO</h2>' + workflow_panel(workflow_from_protocol_graph(g, 'minus_ready')) +
            '<h2>Plus FTO</h2>' + workflow_panel(workflow_from_protocol_graph(g, 'plus_ready')) +
            '<p>Read-through/terminal ACA counting, RNAfold filters and FTO-associated candidate selection '
            'are computational interpretation, not physical RNA removal. The separate ChpBK assays are not modeled as a sequencing branch.</p></main></html>')


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    source = Source(sources_root() / HERE.parent.name / 'sciadv.aax0250_PMC6609220.xml.txt',
                    'Published RNA materials, RNA-ACA control')
    if have(c, source, label='m6A-REF source control'):
        c('validation RNA-ACA is transcribed from its own original', CONTROL_RNA in source.path.read_text(), True)
    treated = RNAHistory().demethylation_treatment().mazf_treatment().repair()
    c('FTO history remains attached through digestion/repair', treated.fto_exposed, True)
    c('matched minus-enzyme route is not silently demethylated', RNAHistory().mazf_treatment().repair().fto_exposed, False)
    for label, operation in [('repair before digest', lambda: RNAHistory().repair()),
                             ('FTO after MazF', lambda: RNAHistory().mazf_treatment().demethylation_treatment()),
                             ('repeat MazF', lambda: RNAHistory().mazf_treatment().mazf_treatment())]:
        try:
            operation()
        except ValueError:
            rejected = True
        else:
            rejected = False
        c('represented route rejects ' + label, rejected, True)
    c('unreviewed historical kit blocks finished publication', protocol_graph().publishable, False)
    c.report()


if __name__ == '__main__':
    selftest()
