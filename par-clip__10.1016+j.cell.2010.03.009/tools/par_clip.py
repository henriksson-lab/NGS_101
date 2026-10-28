"""Original PAR-CLIP front end; no guessed library adapters or mutation chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE = 'Hafner et al. 2010 doi:10.1016/j.cell.2010.03.009, Methods PAR-CLIP (JATS205)'
BOUNDARY = ('Original chemistry supplement is gated. Library preparation explicitly cites '
            'Hafner2008 (10.1016/j.ymeth.2007.09.009), whose Solexa route is read, but '
            'source-specific adapter/variant applicability and downstream molecular model '
            'remain deferred. No sequencing-ready terminal is declared.')


def rna(name):
    return Construct([Segment(name, '[variable 4SU-labelled RNA]', placeholder=True)], name=name)


def protocol_graph():
    g = ProtocolGraph('PAR-CLIP — 4SU crosslink and RNA recovery', authority=SOURCE,
                      representation='published_instance')
    sample = g.input('labelled_cells', 'HEK293 cells grown overnight with 100 µM 4SU',
                     MaterialKind.BIOLOGICAL, primary=True)
    uv = g.input('uv365', '365 nm UV irradiation', MaterialKind.CONDITION)
    crosslinked = g.product('crosslinked_rnp', '4SU-containing RNA crosslinked to RBP',
                            MaterialKind.MOLECULAR, state=rna('crosslinked RNA'))
    g.react('crosslink', 'Irradiate live 4SU-labelled cells at 365 nm', StepKind.CROSS_LINKING,
            ((sample, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (crosslinked,))
    nt1 = g.input('rnase_t1', 'RNase T1', MaterialKind.CATALYST)
    trimmed = g.product('trimmed_rnp', 'RNase-T1-trimmed crosslinked RNA–protein complexes',
                        MaterialKind.MOLECULAR, state=rna('trimmed RNA'))
    g.react('digest_rna', 'NP40 lysis and RNase T1 digestion', StepKind.DIGESTION,
            ((crosslinked, InputRole.TEMPLATE), (nt1, InputRole.CATALYST)), (trimmed,))
    antibody = g.input('anti_flag', 'anti-FLAG Protein G Dynabeads', MaterialKind.SELECTION)
    ip = g.product('ip_rnp', 'immunoprecipitated FLAG/HA-tagged RBP–RNA complexes',
                   MaterialKind.MOLECULAR, state=rna('immunopurified RNA'))
    reject = g.product('ip_waste', 'unretained lysate and wash material', MaterialKind.WASTE)
    g.react('immunoprecipitate', 'Anti-FLAG immunoprecipitation and on-bead RNase T1 digestion',
            StepKind.AFFINITY_CAPTURE, ((trimmed, InputRole.TEMPLATE),
                                      (antibody, InputRole.SELECTION), (nt1, InputRole.CATALYST)), (ip, reject))
    g.discard(reject, 'Immunoprecipitation flow-through and wash fractions')
    cip = g.input('cip', 'calf intestinal alkaline phosphatase', MaterialKind.CATALYST)
    dephos = g.product('dephosphorylated', 'dephosphorylated complex-associated RNA',
                      MaterialKind.MOLECULAR, state=rna('dephosphorylated RNA'))
    g.react('dephosphorylate', 'CIP dephosphorylation on beads', StepKind.DEPHOSPHORYLATION,
            ((ip, InputRole.TEMPLATE), (cip, InputRole.CATALYST)), (dephos,))
    pnk = g.input('pnk', 'polynucleotide kinase', MaterialKind.CATALYST)
    atp = g.input('radioactive_atp', 'radioactive ATP', MaterialKind.CONDITION)
    labelled = g.product('labelled_rnp', 'radiolabelled crosslinked RNA–protein complex',
                         MaterialKind.MOLECULAR, state=rna('radiolabelled RNA'))
    g.react('radiolabel', 'PNK and radioactive ATP label complex-associated RNA', StepKind.PHOSPHORYLATION,
            ((dephos, InputRole.TEMPLATE), (pnk, InputRole.CATALYST), (atp, InputRole.CONDITION)), (labelled,))
    pk = g.input('proteinase_k', 'proteinase K', MaterialKind.CATALYST)
    recovered = g.product('recovered_rna', 'RNA recovered after gel/electroelution and proteinase K extraction',
                           MaterialKind.MOLECULAR, state=rna('recovered RNA'))
    g.react('recover', 'SDS-PAGE, electroelution, proteinase K, acidic phenol/chloroform and precipitation',
            StepKind.DEPROTEINIZATION, ((labelled, InputRole.TEMPLATE), (pk, InputRole.CATALYST)), (recovered,))
    g.mark_terminal(Terminal(recovered.id, Disposition.UNRESOLVED, BOUNDARY))
    g.declare_reviewed_scope(('labelled_cells','uv365','crosslink','rnase_t1','digest_rna','anti_flag',
                             'immunoprecipitate','cip','dephosphorylate','pnk','radioactive_atp',
                             'radiolabel','proteinase_k','recover'), SOURCE)
    return g


def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('recovery endpoint is unresolved, never an invented library',
      g.terminals['recovered_rna'].disposition, Disposition.UNRESOLVED)
    c('source-limited graph fails sequencing publication gate',g.publishable,False)
    c.report()

if __name__ == '__main__': selftest()
