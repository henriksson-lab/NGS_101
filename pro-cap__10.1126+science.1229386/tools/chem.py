from nascent_protocols import make_spec, role_library

from nascent_protocols import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, molecule, role
DEFINITION = Definition("PRO-cap — capped precision run-on sequencing", "10.1126/science.1229386", "Kwak et al., Science (2013)", "Precision run-on with four biotin-NTPs labels nascent RNA while cap selection redirects the library to transcription start sites.", "The cap-selected read begins at the initiating nucleotide; PRO-seq from the same pool reports the active-site 3′ end.", "01_pro_cap.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))


def _ligate_adapter(graph, end, template, state, reviewed):
    adapter = graph.input(end + '_adapter', end + ' RNA adapter', MaterialKind.ADAPTER)
    ligase = graph.input(end + '_ligase', end + ' RNA ligase', MaterialKind.CATALYST)
    product = graph.product(end + '_product', end + '-adapter-ligated nascent RNA',
                            MaterialKind.MOLECULAR, state=state)
    graph.react(end + '_ligation', 'Ligate the ' + end + ' RNA adapter',
                StepKind.LIGATION,
                ((template, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (product,),
                geometry=AdapterLigationGeometry.unavailable(
                    'The 2013 method names standard Illumina RNA adapters but does not print the exact ligating-end sequence used in this experiment.'))
    reviewed.extend((end + '_adapter', end + '_ligase', end + '_ligation'))
    return product


def _runon_cap_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f'{d.citation}; doi:{d.doi}'
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input('nuclei', 'isolated nuclei containing engaged RNA polymerases', MaterialKind.BIOLOGICAL, state=molecule('isolated nuclei', (role('engaged RNA polymerase'), role('nascent RNA'))), primary=True)
    label_name = 'four biotin-NTPs'
    label = graph.input('labelled_ntp', label_name + ' run-on substrates', MaterialKind.MOLECULAR)
    polymerase = graph.input('engaged_polymerase', 'engaged endogenous RNA polymerase', MaterialKind.CATALYST)
    runon = graph.product('runon_rna', 'affinity-labelled nuclear run-on RNA', MaterialKind.MOLECULAR, state=molecule('labelled run-on RNA', (role('capped nascent 5-prime end'), role('labelled run-on 3-prime end'))))
    graph.react('run_on', 'Extend engaged polymerases with ' + label_name, StepKind.RUN_ON, ((nuclei, InputRole.TEMPLATE), (label, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (runon,))
    capture_label = 'streptavidin affinity reagent'
    capture = graph.input('runon_capture', capture_label, MaterialKind.SELECTION)
    enriched = graph.product('enriched_runon', 'affinity-enriched run-on RNA', MaterialKind.MOLECULAR, state=runon.state)
    unlabelled = graph.product('unlabelled_rna', 'RNA not extended in the run-on', MaterialKind.WASTE)
    graph.react('enrich_runon', 'Affinity-enrich labelled run-on RNA', StepKind.AFFINITY_CAPTURE, ((runon, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'did not acquire the run-on affinity label')
    reviewed = ['nuclei', 'labelled_ntp', 'engaged_polymerase', 'run_on', 'runon_capture', 'enrich_runon']
    ligated3 = _ligate_adapter(graph, 'three_prime', enriched,
                              molecule('3-prime-adapter-ligated capped RNA',
                                       (role('capped nascent 5-prime end'), role('nascent RNA body'),
                                        role('3-prime RNA adapter'))), reviewed)
    enriched3 = graph.product('enriched_three_prime_ligated_rna', 'biotin-enriched 3-prime-adapter-ligated RNA', MaterialKind.MOLECULAR, state=ligated3.state)
    rejected3 = graph.product('rejected_three_prime_ligated_rna', 'unlabelled RNA removed after 3-prime ligation', MaterialKind.WASTE)
    graph.react('enrich_three_prime_ligated_rna', 'Re-enrich biotin-labelled RNA after 3-prime ligation', StepKind.AFFINITY_CAPTURE, ((ligated3, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched3, rejected3))
    graph.discard(rejected3, 'lacks the run-on biotin label')
    reviewed.append('enrich_three_prime_ligated_rna')
    phosphatase = graph.input('phosphatase', 'uncapped-RNA 5-prime phosphatase', MaterialKind.CATALYST)
    dephosphorylated = graph.product('dephosphorylated_rna', '3-prime-ligated RNA after uncapped-end blocking', MaterialKind.MOLECULAR, state=enriched3.state)
    graph.react('dephosphorylate', 'Remove 5-prime phosphates from uncapped RNA', StepKind.DEPHOSPHORYLATION, ((enriched3, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,))
    tap = graph.input('tap', 'tobacco acid pyrophosphatase', MaterialKind.CATALYST)
    cap_selected = graph.product('cap_selected_rna', 'decapped nascent RNA with ligatable 5-prime monophosphate', MaterialKind.MOLECULAR, state=molecule('PRO-cap cap-selected RNA', (role('ligatable decapped 5-prime monophosphate'), role('nascent RNA body'), role('3-prime RNA adapter'))))
    graph.react('decap', 'Decap protected RNA ends with TAP', StepKind.CHEMICAL_CONVERSION, ((dephosphorylated, InputRole.TEMPLATE), (tap, InputRole.CATALYST)), (cap_selected,))
    reviewed.extend(['phosphatase', 'dephosphorylate', 'tap', 'decap'])
    ligated5 = _ligate_adapter(graph, 'five_prime', cap_selected,
                              molecule('PRO-cap RNA with both adapters',
                                       (role('5-prime RNA adapter'), role('nascent RNA body'),
                                        role('3-prime RNA adapter'))), reviewed)
    enriched5 = graph.product('enriched_five_prime_ligated_rna', 'biotin-enriched RNA with both adapters', MaterialKind.MOLECULAR, state=ligated5.state)
    rejected5 = graph.product('rejected_five_prime_ligated_rna', 'unlabelled RNA removed after 5-prime ligation', MaterialKind.WASTE)
    graph.react('enrich_five_prime_ligated_rna', 'Re-enrich biotin-labelled RNA after 5-prime ligation', StepKind.AFFINITY_CAPTURE, ((ligated5, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched5, rejected5))
    graph.discard(rejected5, 'lacks the run-on biotin label')
    reviewed.append('enrich_five_prime_ligated_rna')
    finish_cdna(graph, 'library', enriched5, final, d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return _runon_cap_graph(d, final)
