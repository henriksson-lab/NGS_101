from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, ligate_rna_adapter, molecule, role
DEFINITION = Definition("PRO-cap — capped precision run-on sequencing", "10.1126/science.1229386", "Kwak et al., Science (2013)", "Precision run-on marks engaged polymerases with one biotin-NTP while cap selection redirects the library to nascent transcription start sites.", "The cap-selected read begins at the initiating nucleotide; PRO-seq from the same pool reports the active-site 3′ end.", "01_pro_cap.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))

def _runon_cap_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f'{d.citation}; doi:{d.doi}'
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input('nuclei', 'isolated nuclei containing engaged RNA polymerases', MaterialKind.BIOLOGICAL, state=molecule('isolated nuclei', (role('engaged RNA polymerase'), role('nascent RNA'))), primary=True)
    label_name = 'single biotin-NTP'
    label = graph.input('labelled_ntp', label_name + ' run-on nucleotide', MaterialKind.MOLECULAR)
    polymerase = graph.input('engaged_polymerase', 'engaged endogenous RNA polymerase', MaterialKind.CATALYST)
    runon = graph.product('runon_rna', 'affinity-labelled nuclear run-on RNA', MaterialKind.MOLECULAR, state=molecule('labelled run-on RNA', (role('capped nascent 5-prime end'), role('labelled run-on 3-prime end'))))
    graph.react('run_on', 'Extend engaged polymerases with ' + label_name, StepKind.RUN_ON, ((nuclei, InputRole.TEMPLATE), (label, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (runon,))
    capture_label = 'streptavidin affinity reagent'
    capture = graph.input('runon_capture', capture_label, MaterialKind.SELECTION)
    enriched = graph.product('enriched_runon', 'affinity-enriched run-on RNA', MaterialKind.MOLECULAR, state=runon.state)
    unlabelled = graph.product('unlabelled_rna', 'RNA not extended in the run-on', MaterialKind.WASTE)
    graph.react('enrich_runon', 'Affinity-enrich labelled run-on RNA', StepKind.AFFINITY_CAPTURE, ((runon, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'did not acquire the run-on affinity label')
    phosphatase = graph.input('phosphatase', 'uncapped-RNA 5-prime phosphatase', MaterialKind.CATALYST)
    dephosphorylated = graph.product('dephosphorylated_rna', 'run-on RNA after uncapped-end blocking', MaterialKind.MOLECULAR, state=enriched.state)
    graph.react('dephosphorylate', 'Remove 5-prime phosphates from uncapped RNA', StepKind.DEPHOSPHORYLATION, ((enriched, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    tap = graph.input('tap', 'tobacco acid pyrophosphatase', MaterialKind.CATALYST)
    cap_selected = graph.product('cap_selected_rna', 'decapped nascent RNA with ligatable 5-prime monophosphate', MaterialKind.MOLECULAR, state=molecule('PRO-cap' + ' cap-selected RNA', (role('ligatable decapped 5-prime monophosphate'), role('nascent RNA body'), role('labelled run-on 3-prime end'))))
    graph.react('decap', 'Decap protected RNA ends with TAP', StepKind.CHEMICAL_CONVERSION, ((dephosphorylated, InputRole.TEMPLATE), (tap, InputRole.CATALYST)), (cap_selected,))
    reviewed = ['nuclei', 'labelled_ntp', 'engaged_polymerase', 'run_on', 'runon_capture', 'enrich_runon', 'phosphatase', 'dephosphorylate', 'tap', 'decap']
    ligated5 = ligate_rna_adapter(graph, 'five_prime', cap_selected, '5-prime', reviewed)
    ligated3 = ligate_rna_adapter(graph, 'three_prime', ligated5, '3-prime', reviewed)
    finish_cdna(graph, 'library', ligated3, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return _runon_cap_graph(d, final)
