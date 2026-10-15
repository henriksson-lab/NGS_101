from nascent_protocols import make_spec, role_library

from adapter_ligation import reported_adapter_end
from chemdraw import Segment
from illumina import TRUSEQ_SMALL_RNA_RA3, TRUSEQ_SMALL_RNA_RA5
from nascent_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, molecule, role
from protocol_model import molecule as drawn_molecule

KRUESI = 'Kruesi et al., eLife 2013, doi:10.7554/eLife.00808, Methods pp. 22–23'
ILLUMINA = 'Illumina TruSeq Small RNA adapter sequences, document 1000000002694 v22'
RA3 = TRUSEQ_SMALL_RNA_RA3
RA5 = TRUSEQ_SMALL_RNA_RA5
DEFINITION = Definition("GRO-cap — capped global run-on sequencing", "10.1038/ng.3142", "Core et al., Nature Genetics (2014)", "BrUTP-labelled global run-on RNA is subjected to enzymatic cap selection, enriching the capped 5′ ends of nascent promoter and enhancer transcripts.", "The informative base is the first transcribed nucleotide immediately downstream of the ligated 5′ adapter.", "01_gro_cap.html", "Kruesi et al. (2013), explicitly followed by Core et al. (2014), ligate RA3 before cap selection and RA5 afterward. The published no-TAP aliquot measures background; this schematic follows the TAP+ library. Illumina's current guide prints the named RA3 and RA5 bases, but equivalence to the 2013 kit revision and the exact selected index remain unverified.")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The selected six-base index and exact historical primer revision are unresolved; Core et al. 2008 supplementary methods identify the inherited RT and PCR operations.',
    oligos=lambda: (
        f'<b>RA3, named by Kruesi; current Illumina kit reference</b> — 5′-<inf>{RA3}</inf>-3′ (RNA 3′ adapter; vendor uses T spelling)',
        f'<b>RA5, named by Kruesi; current Illumina kit reference</b> — 5′-<inf>{RA5.replace("T", "U")}</inf>-3′ (vendor RNA spelling; model uses T for construct arithmetic)',
    ))))


def _ligate_named_adapter(graph, prefix, template, position, sequence, reviewed):
    """Keep source-named adapter identity and ligating end beside this protocol."""
    label = f'TruSeq RNA {position} adapter ({prefix.upper()})'
    adapter = graph.input(prefix + '_adapter', label, MaterialKind.ADAPTER,
                          state=drawn_molecule(label, (Segment(label, sequence, inferred=True),)))
    ligase = graph.input(prefix + '_ligase', 'T4 RNA ligase I', MaterialKind.CATALYST)
    end = '5-prime' if position == '3-prime' else '3-prime'
    geometry = reported_adapter_end(
        label, kind='single-strand', polarity=end, terminal_sequence=sequence,
        statement=f'{KRUESI} names {label} and its ligation position; {ILLUMINA} prints this sequence. Historical kit revision is not independently confirmed.')
    if position == '3-prime':
        parts = (role('protected capped 5-prime end'), role('nascent RNA body'),
                 Segment('RA3', sequence, inferred=True))
        junctions = (('nascent RNA body', 'RA3', 'ligation'),)
    else:
        parts = (Segment('RA5', sequence, inferred=True), role('nascent RNA body'),
                 Segment('RA3', RA3, inferred=True))
        junctions = (('RA5', 'nascent RNA body', 'ligation'),
                     ('nascent RNA body', 'RA3', 'ligation'))
    product = graph.product(prefix + '_product', label + '-ligated RNA',
                            MaterialKind.MOLECULAR,
                            state=drawn_molecule(label + '-ligated RNA', parts,
                                                 junctions=junctions))
    graph.react(prefix + '_ligation', 'Ligate ' + label,
                StepKind.LIGATION,
                ((template, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (product,),
                geometry=geometry, source=KRUESI)
    reviewed.extend((prefix + '_adapter', prefix + '_ligase', prefix + '_ligation'))
    return product

def _runon_cap_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f'{d.citation}; doi:{d.doi}'
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input('nuclei', 'isolated nuclei containing engaged RNA polymerases', MaterialKind.BIOLOGICAL, state=molecule('isolated nuclei', (role('engaged RNA polymerase'), role('nascent RNA'))), primary=True)
    label_name = 'BrUTP'
    label = graph.input('labelled_ntp', label_name + ' run-on nucleotide', MaterialKind.MOLECULAR)
    polymerase = graph.input('engaged_polymerase', 'engaged endogenous RNA polymerase', MaterialKind.CATALYST)
    runon = graph.product('runon_rna', 'affinity-labelled nuclear run-on RNA', MaterialKind.MOLECULAR, state=molecule('labelled run-on RNA', (role('capped nascent 5-prime end'), role('labelled run-on 3-prime end'))))
    graph.react('run_on', 'Extend engaged polymerases with ' + label_name, StepKind.RUN_ON, ((nuclei, InputRole.TEMPLATE), (label, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (runon,))
    capture_label = 'anti-BrdU affinity reagent'
    capture = graph.input('runon_capture', capture_label, MaterialKind.SELECTION)
    enriched = graph.product('enriched_runon', 'affinity-enriched run-on RNA', MaterialKind.MOLECULAR, state=runon.state)
    unlabelled = graph.product('unlabelled_rna', 'RNA not extended in the run-on', MaterialKind.WASTE)
    graph.react('enrich_runon', 'Affinity-enrich labelled run-on RNA', StepKind.AFFINITY_CAPTURE, ((runon, InputRole.TEMPLATE), (capture, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'did not acquire the run-on affinity label')
    reviewed = ['nuclei', 'labelled_ntp', 'engaged_polymerase', 'run_on', 'runon_capture', 'enrich_runon']
    ligated3 = _ligate_named_adapter(graph, 'ra3', enriched, '3-prime', RA3, reviewed)
    capture2 = graph.input('second_capture', 'second anti-BrdU affinity capture', MaterialKind.SELECTION)
    enriched2 = graph.product('twice_enriched_rna', 'RA3-ligated run-on RNA after second BrU capture', MaterialKind.MOLECULAR, state=ligated3.state)
    rejected2 = graph.product('second_capture_reject', 'RNA not retained by second BrU capture', MaterialKind.WASTE)
    graph.react('enrich_after_ra3', 'Affinity-enrich RA3-ligated run-on RNA', StepKind.AFFINITY_CAPTURE,
                ((ligated3, InputRole.TEMPLATE), (capture2, InputRole.SELECTION)),
                (enriched2, rejected2), source=KRUESI, changes_molecule=False)
    graph.discard(rejected2, 'not recovered by the second BrU affinity capture')
    terminator = graph.input('terminator', 'Terminator 5-prime phosphate-dependent exonuclease', MaterialKind.CATALYST)
    terminator_selected = graph.product('terminator_selected', 'capped RNA retained after 5-prime-monophosphate digestion', MaterialKind.MOLECULAR, state=enriched2.state)
    degraded_uncapped = graph.product('degraded_uncapped', 'uncapped monophosphate RNA degraded by Terminator', MaterialKind.WASTE)
    graph.react('terminator_digest', 'Degrade uncapped 5-prime-monophosphate RNA', StepKind.DIGESTION,
                ((enriched2, InputRole.TEMPLATE), (terminator, InputRole.CATALYST)),
                (terminator_selected, degraded_uncapped), source=KRUESI)
    graph.discard(degraded_uncapped, 'Terminator exonuclease removes uncapped 5-prime-monophosphate RNA')
    phosphatase = graph.input('phosphatase', 'Antarctic phosphatase', MaterialKind.CATALYST)
    dephosphorylated = graph.product('dephosphorylated_rna', 'RNA after residual uncapped-end dephosphorylation', MaterialKind.MOLECULAR, state=terminator_selected.state)
    graph.react('dephosphorylate', 'Dephosphorylate residual uncapped 5-prime ends', StepKind.DEPHOSPHORYLATION, ((terminator_selected, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)), (dephosphorylated,), source=KRUESI)
    tap = graph.input('tap', 'tobacco acid pyrophosphatase', MaterialKind.CATALYST)
    cap_selected = graph.product('cap_selected_rna', 'decapped nascent RNA with ligatable 5-prime monophosphate', MaterialKind.MOLECULAR, state=drawn_molecule('GRO-cap cap-selected RNA', (role('ligatable decapped 5-prime monophosphate'), role('nascent RNA body'), Segment('RA3', RA3, inferred=True)), junctions=(('nascent RNA body', 'RA3', 'ligation'),)))
    graph.react('decap', 'Decap protected RNA ends with TAP', StepKind.CHEMICAL_CONVERSION, ((dephosphorylated, InputRole.TEMPLATE), (tap, InputRole.CATALYST)), (cap_selected,))
    reviewed.extend(('second_capture', 'enrich_after_ra3', 'terminator',
                     'terminator_digest', 'phosphatase', 'dephosphorylate',
                     'tap', 'decap'))
    ligated5 = _ligate_named_adapter(graph, 'ra5', cap_selected, '5-prime', RA5, reviewed)
    capture3 = graph.input('third_capture', 'third anti-BrdU affinity capture', MaterialKind.SELECTION)
    enriched3 = graph.product('fully_ligated_enriched', 'RA5/RA3-ligated RNA after third BrU capture', MaterialKind.MOLECULAR, state=ligated5.state)
    rejected3 = graph.product('third_capture_reject', 'RNA not retained by third BrU capture', MaterialKind.WASTE)
    graph.react('enrich_after_ra5', 'Affinity-enrich RA5/RA3-ligated run-on RNA', StepKind.AFFINITY_CAPTURE,
                ((ligated5, InputRole.TEMPLATE), (capture3, InputRole.SELECTION)),
                (enriched3, rejected3), source=KRUESI, changes_molecule=False)
    graph.discard(rejected3, 'not recovered by the third BrU affinity capture')
    reviewed.extend(('third_capture', 'enrich_after_ra5'))
    finish_cdna(graph, 'library', enriched3, final, d.endpoint, reviewed,
                rt_label='TruSeq Small RNA RTP (Core 2008 inherited operation; historical primer revision unresolved)')
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return _runon_cap_graph(d, final)
