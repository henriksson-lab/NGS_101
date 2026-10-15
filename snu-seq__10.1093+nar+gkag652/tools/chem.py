from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, FourSUReadout, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, ligate_rna_adapter, molecule, role
DEFINITION = Definition("SNU-Seq — single-nucleotide-resolution 4sU sequencing", "10.1093/nar/gkag652", "Gerlevik et al., Nucleic Acids Research (2026)", "Biotin-selected intact 4sU RNA is polyadenylated with bacterial poly(A) polymerase and captured by a 3′ RNA-seq workflow to preserve native endpoint information.", "The 3′-anchored read reports the native end of an intact pulse-labelled RNA molecule.", "01_snu_seq.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=False, platform="platform-specific")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))

READOUT = FourSUReadout('thiol biotinylation, streptavidin purification and bacterial poly(A) polymerase', 'native RNA ends', True)

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
    biotin = graph.input('thiol_biotin', '4sU thiol-biotinylation reagent', MaterialKind.MOLECULAR)
    condition = graph.input('biotinylation_condition', 'thiol-biotinylation condition', MaterialKind.CONDITION)
    biotinylated = graph.product('biotinylated_rna', 'biotinylated 4sU RNA', MaterialKind.MOLECULAR, state=labelled.state)
    graph.react('biotinylate', 'Thiol-biotinylate incorporated 4sU', StepKind.CHEMICAL_CONVERSION, ((labelled, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL), (condition, InputRole.CONDITION)), (biotinylated,))
    streptavidin = graph.input('streptavidin', 'streptavidin RNA-capture reagent', MaterialKind.SELECTION)
    enriched = graph.product('enriched_rna', 'streptavidin-enriched 4sU RNA', MaterialKind.MOLECULAR, state=cfg.product())
    unlabelled = graph.product('unlabelled_rna', 'unlabelled RNA background', MaterialKind.WASTE)
    graph.react('enrich', 'Affinity-purify biotinylated 4sU RNA', StepKind.AFFINITY_CAPTURE, ((biotinylated, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'does not carry the pulse-installed 4sU affinity label')
    reviewed.extend(('thiol_biotin', 'biotinylation_condition', 'biotinylate', 'streptavidin', 'enrich'))
    pap = graph.input('poly_a_polymerase', 'bacterial poly(A) polymerase', MaterialKind.CATALYST)
    atp = graph.input('tailing_atp', 'ATP for RNA polyadenylation', MaterialKind.MOLECULAR)
    selected = graph.product('selected_rna', 'polyadenylated intact 4sU RNA', MaterialKind.MOLECULAR, state=enriched.state)
    graph.react('polyadenylate', 'Polyadenylate intact enriched RNA', StepKind.TAILING, ((enriched, InputRole.TEMPLATE), (atp, InputRole.STRUCTURAL), (pap, InputRole.CATALYST)), (selected,))
    reviewed.extend(('poly_a_polymerase', 'tailing_atp', 'polyadenylate'))
    ligated = ligate_rna_adapter(graph, 'library_adapter', selected, 'directional library', reviewed)
    finish_cdna(graph, 'library', ligated, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return four_su_graph(d, final)
