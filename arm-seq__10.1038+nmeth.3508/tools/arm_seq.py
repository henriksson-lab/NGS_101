"""ARM-seq pretreatment graph; historical kit library remains an explicit boundary."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            Terminal, Disposition)
SOURCE = 'Cozen et al. 2015 doi:10.1038/nmeth.3508, Methods Treatment of RNA with AlkB; library preparation'
BOUNDARY = ('The source names NEBNext Small RNA Library Prep Kit for Illumina; applicable '
            'historical kit identity/revision and oligos have not been reviewed. '
            'MiSeq paired-end output is reported, but no exact library or run lengths are derived.')


def rna_state(name):
    # A variable population, not a fabricated example tRNA or residue-specific modification map.
    return Construct([Segment(name, '[variable RNA population]', placeholder=True)], name=name)


def protocol_graph():
    g = ProtocolGraph('ARM-seq — AlkB pretreatment', authority=SOURCE,
                      representation='published_instance')
    rna = g.input('bulk_rna', 'bulk RNA containing AlkB-sensitive modifications',
                  MaterialKind.MOLECULAR, primary=True, state=rna_state('bulk RNA'))
    alkb = g.input('alkb', 'wild-type E. coli AlkB', MaterialKind.CATALYST)
    conditions = g.input('reaction_buffer', 'HEPES-KOH, Fe(II), alpha-ketoglutarate, ascorbate and BSA',
                         MaterialKind.CONDITION)
    treated = g.product('treated_rna', 'AlkB-treated RNA; modification-sensitive RT blocks reduced',
                        MaterialKind.MOLECULAR, state=rna_state('AlkB-treated RNA'))
    g.react('demethylate', 'AlkB demethylation (37 °C, 100 min; source-specific buffer)',
            StepKind.CHEMICAL_CONVERSION, ((rna, InputRole.TEMPLATE), (alkb, InputRole.CATALYST),
                                         (conditions, InputRole.CONDITION)), (treated,))
    clean = g.product('recovered_rna', 'EDTA-quenched, extracted and ethanol-precipitated RNA',
                      MaterialKind.MOLECULAR, state=rna_state('recovered RNA'))
    waste = g.product('cleanup_waste', 'discarded extraction phases and precipitation supernatant', MaterialKind.WASTE)
    g.react('recover', 'EDTA/ammonium-acetate stop, phenol extraction and ethanol precipitation',
            StepKind.CLEANUP, ((treated, InputRole.TEMPLATE),), (clean, waste))
    g.discard(waste, 'Physical extraction and precipitation waste')
    small = g.product('small_rna', 'MirVana-selected <200 nt RNA after DNase I and cleanup',
                      MaterialKind.MOLECULAR, state=rna_state('small RNA'))
    selection = g.input('mirvana', 'MirVana miRNA Isolation Kit', MaterialKind.SELECTION)
    rejected = g.product('size_reject', 'material rejected during small-RNA selection and cleanup', MaterialKind.WASTE)
    g.react('prepare_small_rna', 'MirVana selection, concentration, DNase I and column cleanup',
            StepKind.SELECTION, ((clean, InputRole.TEMPLATE), (selection, InputRole.SELECTION)), (small, rejected))
    g.discard(rejected, 'Physical size-selection and cleanup reject')
    g.mark_terminal(Terminal(small.id, Disposition.UNRESOLVED, BOUNDARY))
    g.declare_reviewed_scope(('bulk_rna', 'alkb', 'reaction_buffer', 'demethylate',
                             'recover', 'prepare_small_rna'), SOURCE)
    return g


def selftest():
    from checks import Check
    c = Check()
    g = protocol_graph()
    c('unreviewed kit does not become a finished library', g.publishable, False)
    c('kit boundary retains molecular RNA endpoint',
      g.terminals['small_rna'].disposition, Disposition.UNRESOLVED)
    c.report()

if __name__ == '__main__':
    selftest()
