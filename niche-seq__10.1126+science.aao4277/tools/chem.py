"""Source-qualified NICHE-seq naive lymph-node frontend; private partial graph.

The inspected original is an associated author manuscript. Its directly cited
MARS-seq component is opaque here; no predecessor chemistry is imported.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, SequencingPlan, SourceItem,
                            Terminal, Disposition)
PHOTOACTIVATION_NM = 820
IMAGING_NM = 940
AUTHORITY = ('Associated author manuscript PMC7234837/EMS86422, main47/49/51, '
             'reference7 and reference17; detailed supplementary Methods gated. '
             'Not asserted equivalent to final published Methods.')

def protocol_graph():
    g = ProtocolGraph('NICHE-seq: naive inguinal lymph-node niche selection',
                      authority=AUTHORITY, representation='published_instance')
    tissue = g.input('prepared_ln', 'Naive inguinal lymph node in prepared PA-GFP host mouse, '
                     'with transferred TdTomato+ T-cell and CFP+ B-cell landmarks; '
                     'upstream transfers/imaging preparation outside selected boundary',
                     K.BIOLOGICAL, primary=True)
    laser = g.input('activation_condition', f'Source two-photon irradiation at {PHOTOACTIVATION_NM}nm '
                    'of selected B-cell follicle or T-cell region; duration/dose unassembled', K.CONDITION)
    activated = g.product('activated_tissue', 'Source spatially activated PA-GFP protein '
                          'reporter population in tissue; no nucleotide barcode or perfect ROI confinement asserted', K.MOLECULAR)
    g.react('photoactivate', 'Source local protein-reporter photoactivation', S.CHEMICAL_CONVERSION,
            ((tissue, R.TEMPLATE), (laser, R.CONDITION)), (activated,))
    imaging = g.input('imaging_condition', f'Source subsequent imaging at {IMAGING_NM}nm; '
                      'optical readout of reporter, no sequencing cycles', K.CONDITION)
    imaged = g.product('imaged_tissue', 'Same prepared tissue after optical inspection', K.MOLECULAR)
    g.react('image', 'Source confirm activated fluorescence', S.GENERIC,
            ((activated, R.TEMPLATE), (imaging, R.CONDITION)), (imaged,), changes_molecule=False)
    dissociation = g.input('dissociation_condition', 'Source tissue dissociation; exact historical '
                           'enzymes, doses, fractions and viability gate unassembled from unread supplementary Methods', K.CONDITION)
    cells = g.product('cell_suspension', 'Source dissociated tissue cell suspension; '
                      'recovery and remaining tissue fractions unassembled', K.MOLECULAR)
    g.react('dissociate', 'Opaque source tissue dissociation boundary', S.MAGIC,
            ((imaged, R.TEMPLATE), (dissociation, R.CONDITION)), (cells,),
            reason='Main47 confirms dissociation before FACS, but detailed supplementary Methods are gated; no invented recovery or physical fractions.')
    gate = g.input('facs_condition', 'Source fluorescence-based FACS of photoactivated PA-GFP+ '
                   'single cells; full gates/plate stocks unassembled', K.CONDITION)
    selected = g.product('sorted_cells', 'Physically sorted photoactivated single-cell preparation; '
                         'not a computational subset of total-tissue reads', K.MOLECULAR)
    excluded = g.product('other_cells', 'Unselected cells outside chosen FACS collection; '
                         'amount and later fate unassigned', K.MOLECULAR)
    g.react('facs', 'Source physical single-cell FACS collection', S.SAMPLING,
            ((cells, R.TEMPLATE), (gate, R.CONDITION)), (selected, excluded))
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED,
                            'Source selection does not establish the later fate of uncollected cells.'))
    stocks = g.input('mars_stocks', 'Source-delegated MARS-seq library-preparation stocks; '
                     'exact local historical version, quantities and whole oligos unassembled', K.MOLECULAR)
    library = g.product('library', 'Source MARS-seq library from physically selected niche cells; '
                        'complete molecular states unassembled', K.LIBRARY)
    g.react('mars', 'Source MARS-seq sequencing preparation, opaque delegated component', S.MAGIC,
            ((selected, R.TEMPLATE), (stocks, R.STRUCTURAL)), (library,),
            reason='Main47/51 delegates MARS-seq and reference7 identifies Jaitin2014. Neither its original nor detailed local Methods is inspected; no constructor inheritance or current kit substitution.')
    g.sequence(library, SequencingPlan(limitation='Exact sequencer/read/index masks, cell/UMI barcodes and '
               'primer landing unassembled. PA-GFP is a protein fluorescence label; no inline spatial DNA index '
               'is inferred. Unactivated controls and other tissue/application cohorts are separate comparisons, '
               'not extra collected protocols.'))
    g.declare_reviewed_scope(('prepared_ln','activation_condition','photoactivate','imaging_condition','image','dissociation_condition','dissociate','facs_condition','facs'),AUTHORITY)
    g.add_source_item(SourceItem('pending', AUTHORITY,
                      'Selected frontend only; full supplementary Methods, whole library and final-version applicability remain open. No confirmed absence.'))
    return g

def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    g = protocol_graph()
    c('FACS collection uses dissociated tissue', g.steps['facs'].inputs[0].material, 'cell_suspension')
    c('Delegated preparation consumes physical selected cells', g.steps['mars'].inputs[0].material, 'sorted_cells')
    c('Unselected cells are not silently discarded reads', g.terminals['other_cells'].disposition, Disposition.UNRESOLVED)
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
                    'science.aao4277_PMC7234837.xml.txt', 'https://doi.org/10.1126/science.aao4277')
    if have(c, source):
        text = source.path.read_text()
        c('Whole source activation wavelength', f'at {PHOTOACTIVATION_NM} nm' in text, True)
        c('Whole source imaging wavelength', f'at {IMAGING_NM} nm' in text, True)
        c('Inspected original is author manuscript', 'Funders Author Manuscripts' in text, True)
    c('Partial frontend cannot claim full sequencing readiness', bool(g.publication_problems()), True)
    c.report()

if __name__ == '__main__':
    selftest()
