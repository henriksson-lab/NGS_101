from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, FourSUReadout, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, ligate_rna_adapter, molecule, role
DEFINITION = Definition("TimeLapse-seq — nucleoside-recoding RNA dynamics", "10.1038/nmeth.4582", "Schofield et al., Nature Methods (2018)", "Oxidative nucleophilic substitution converts 4sU in newly synthesized RNA into a cytidine-like base, producing T-to-C mismatches in ordinary RNA-seq reads.", "Standard RNA-seq reads identify newly synthesized molecules from distributed T-to-C conversions.", "01_timelapse_seq.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))

READOUT = FourSUReadout('oxidation plus 2,2,2-trifluoroethylamine nucleophilic substitution', 'U-to-C mismatch', False)

def four_su_graph(d: Definition, final: Construct) -> ProtocolGraph:
    cfg = READOUT
    source = f'{d.citation}; doi:{d.doi}'
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input('cells', 'living cells during the 4sU pulse', MaterialKind.BIOLOGICAL, state=molecule('living cells', (role('RNA synthesis during pulse'),)), primary=True)
    four_su = graph.input('four_su', '4-thiouridine pulse-label nucleotide', MaterialKind.MOLECULAR)
    polymerase = graph.input('cellular_rna_polymerase', 'cellular RNA polymerases', MaterialKind.CATALYST)
    labelled = graph.product('labelled_rna', '4sU pulse-labelled RNA', MaterialKind.MOLECULAR, state=cfg.labelled())
    graph.react('pulse_label', 'Pulse-label newly synthesized RNA with 4sU', StepKind.TRANSCRIPTION, ((cells, InputRole.TEMPLATE), (four_su, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (labelled,))
    reviewed = ['cells', 'four_su', 'cellular_rna_polymerase', 'pulse_label']
    chemistry = graph.input('recoding_chemistry', cfg.chemistry, MaterialKind.CONDITION)
    selected = graph.product('recoded_rna', 'chemically recoded 4sU RNA', MaterialKind.MOLECULAR, state=cfg.product())
    graph.react('recode', 'Chemically recode 4sU before reverse transcription', StepKind.CHEMICAL_CONVERSION, ((labelled, InputRole.TEMPLATE), (chemistry, InputRole.CONDITION)), (selected,))
    reviewed.extend(('recoding_chemistry', 'recode'))
    ligated = ligate_rna_adapter(graph, 'library_adapter', selected, 'directional library', reviewed)
    finish_cdna(graph, 'library', ligated, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return four_su_graph(d, final)
