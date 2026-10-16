from dataclasses import dataclass
from types import SimpleNamespace

from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Segment
from illumina import TRUSEQ_SMALL_RNA_RA3, TRUSEQ_SMALL_RNA_RA5
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from protocol_model import molecule as drawn_molecule

@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    endpoint: str
    notes: str
    caveat: str = ""


def make_spec(d: Definition, final: Construct, *, unavailable: str,
              oligos=lambda: ()):
    """Assemble page fields from choices declared by the owning protocol."""
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint, SEQUENCING_UNAVAILABLE=unavailable,
        oligos=oligos)


KRUESI = 'Kruesi et al., eLife 2013, doi:10.7554/eLife.00808, Methods pp. 22–23'
ILLUMINA = 'Illumina TruSeq Small RNA adapter sequences, document 1000000002694 v22'
RA3 = TRUSEQ_SMALL_RNA_RA3
RA5 = TRUSEQ_SMALL_RNA_RA5
RT_PROGRAM = ((44, 15), (52, 45))  # Core 2008 original SOM p. 19
PCR_CYCLES = 15
PAGE_MIN_NT = 90  # strictly greater than this threshold
DEFINITION = Definition("GRO-cap — capped global run-on sequencing", "10.1038/ng.3142", "Core et al., Nature Genetics (2014)", "BrUTP-labelled global run-on RNA is subjected to enzymatic cap selection, enriching the capped 5′ ends of nascent promoter and enhancer transcripts.", "The informative base is the first transcribed nucleotide immediately downstream of the ligated 5′ adapter.", "01_gro_cap.html", "Kruesi et al. (2013), explicitly followed by Core et al. (2014), ligate RA3 before cap selection and RA5 afterward. The published no-TAP aliquot measures background; this schematic follows the TAP+ library. Illumina's 2011 guide corroborates Kruesi's RA3 and RA5 part numbers, while its current sequence table prints their bases. Kruesi also prints the RA3 read-through trimming target; this is not a full modified ordering strand. Exact historical reagent-sequence continuity and the 2014 GRO-cap primer/index assignment remain unverified; the adjacent GRO-seq RPI family is not assigned here.")

# The GRO-seq paragraph next to GRO-cap prints an RPI six-base-index family.
# That adjacent method does not assign a primer/index to the GRO-cap library.
FINAL_LIBRARY = Construct([
    role("unassembled library left end"), role("GRO-cap assay-derived cDNA"),
    role("unassembled library right end"),
], name=DEFINITION.title + " unresolved physical endpoint")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The historical GRO-cap primer pair, indexing assignment and completed molecular endpoint are not assembled. The adjacent GRO-seq RPI family does not establish a GRO-cap index or its length. Core et al. 2008 identifies the inherited RT, RNA digestion, PCR and PAGE operations.',
    oligos=lambda: (
        f'<b>RA3, named by Kruesi; current Illumina kit reference</b> — 5′-<inf>{RA3}</inf>-3′ (RNA 3′ adapter; vendor uses T spelling)',
        f'<b>RA5, named by Kruesi; current Illumina kit reference</b> — 5′-<inf>{RA5.replace("T", "U")}</inf>-3′ (vendor RNA spelling; model uses T for construct arithmetic)',
    ))))
ACQUIRED_ROLES=("Read 1",)
READ_SPECS=(sp.read_spec("Read 1", "acquired cycles unassigned",
    "Deposited GRO-cap experiments declare single-end HiSeq 2000 acquisition; mapped reads are trimmed to 30 bases.",
    "include"),)
SEQUENCING_ENDING="GEO-linked SRA experiments SRX682016–SRX682019 declare single-end Illumina HiSeq 2000 acquisition. The 30-base mapped/trimmed length is not an acquired cycle count."


def _ligate_named_adapter(graph, prefix, template, position, sequence, reviewed):
    """Keep source-named adapter identity and ligating end beside this protocol."""
    label = f'TruSeq RNA {position} adapter ({prefix.upper()})'
    adapter = graph.input(prefix + '_adapter', label, MaterialKind.ADAPTER,
                          state=drawn_molecule(label, (Segment(label, sequence, inferred=True),)))
    ligase = graph.input(prefix + '_ligase', 'T4 RNA ligase I', MaterialKind.CATALYST)
    end = '5-prime' if position == '3-prime' else '3-prime'
    geometry = reported_adapter_end(
        label, kind='single-strand', polarity=end, terminal_sequence=sequence,
        statement=f'{KRUESI} names {label} and its ligation position; Illumina Rev C (2011) corroborates the component part number; {ILLUMINA} prints this sequence. Exact 2013 reagent-sequence continuity is unverified.')
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
    graph.react(prefix + '_ligation', 'INFERRED — ligate ' + label + ' using current kit-reference bases',
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
    completion_source = 'Core et al. 2008 original Science SOM pp. 18–19, inherited explicitly by Kruesi GRO-cap'
    rt_primer = graph.input('library_rt_primer', 'Illumina Small RNA kit RT primer; historical strand unassigned', MaterialKind.OLIGO, source=completion_source)
    rt_enzyme = graph.input('library_reverse_transcriptase', 'SuperScript III', MaterialKind.CATALYST, source=completion_source)
    hybrid = graph.product('library_cdna', 'reverse-transcribed RNA/cDNA material', MaterialKind.MOLECULAR,
                           state=molecule('GRO-cap RT product', (role('RNA template'), role('first-strand cDNA'))))
    programme = ', then '.join(f'{temperature}°C for {minutes} min' for temperature, minutes in RT_PROGRAM)
    graph.react('library_rt', 'Reverse-transcribe with the inherited kit primer: ' + programme, StepKind.REVERSE_TRANSCRIPTION,
                ((enriched3, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt_enzyme, InputRole.CATALYST)),
                (hybrid,), source=completion_source)
    rnases = graph.input('library_rnases', 'RNase cocktail and RNase H', MaterialKind.CATALYST, source=completion_source)
    cdna = graph.product('library_rna_free_cdna', 'cDNA after RNA-template digestion', MaterialKind.MOLECULAR,
                         state=molecule('RNA-free cDNA', (role('first-strand cDNA'),)))
    digested_rna = graph.product('library_digested_rna', 'digested RNA templates', MaterialKind.WASTE)
    graph.react('library_digest_rna', 'Digest the RNA template before PCR', StepKind.DIGESTION,
                ((hybrid, InputRole.TEMPLATE), (rnases, InputRole.CATALYST)), (cdna, digested_rna), source=completion_source)
    graph.discard(digested_rna, 'RNA digestion products are not the amplified cDNA library')
    forward = graph.input('library_pcr_forward', 'Illumina kit PCR forward primer; historical strand unassigned', MaterialKind.OLIGO, source=completion_source)
    reverse = graph.input('library_pcr_reverse', 'Illumina kit PCR reverse primer; index assignment unresolved', MaterialKind.OLIGO, source=completion_source)
    phusion = graph.input('library_pcr_polymerase', 'Phusion high-fidelity DNA polymerase', MaterialKind.CATALYST, source=completion_source)
    amplified = graph.product('library_amplified', 'PCR-amplified cDNA before PAGE', MaterialKind.MOLECULAR, state=final)
    graph.react('library_pcr', f'Amplify {PCR_CYCLES} cycles with the inherited kit PCR primers', StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)),
                (amplified,), source=completion_source, geometry=PCRGeometry.unavailable(
                    cdna.state, reason='The historical GRO-cap primer pair and its precursor binding sites are not assembled in this model.'))
    selection = graph.input('library_page_selection', f'inherited non-denaturing PAGE selection above {PAGE_MIN_NT} nt', MaterialKind.SELECTION, source=completion_source)
    library = graph.product('library_library', 'PAGE-purified GRO-cap library', MaterialKind.LIBRARY, state=final)
    short = graph.product('library_short_products', 'products outside the selected PAGE range', MaterialKind.WASTE)
    graph.react('library_page', 'PAGE-purify the amplified cDNA', StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (selection, InputRole.SELECTION)), (library, short),
                source=completion_source, changes_molecule=False)
    graph.discard(short, 'products outside the inherited PAGE selection are not sequenced')
    graph.sequence(library, SequencingPlan(library=final, primers=(), acquired_roles=ACQUIRED_ROLES,
                                          limitation=SEQUENCING_UNAVAILABLE))
    graph.annotate(SourceDiscrepancy(
        target="ra3_ligation",title="Deposited protocol volume units differ from the original paper",
        claims=(
            SourceClaim("Kruesi 2013 original PDF Methods pp.22–23, RA3 ligation paragraph",
                        "Prints 2 microlitres of 50% PEG, bringing the mixture to 14 microlitres."),
            SourceClaim("GEO GSE60456/GSE43087 GRO-cap sample extraction protocols and Core2014 ENA construction metadata",
                        "Render the corresponding quantities as 2 ml and 14 ml.")),
        model_choice="Use the original paper for method interpretation; this model assigns neither of those volumes.",
        evidence="The deposited paragraphs otherwise reproduce the same named adapters, reagents and reaction order. The original PDF visibly uses the microlitre symbol.",
        uncertainty="This is an apparent metadata transcription/unit error; no formal source correction is established.",
        slot=AnnotationSlot.INTRO))
    reviewed.extend(('library_rt_primer', 'library_reverse_transcriptase', 'library_rt', 'library_rnases',
                     'library_digest_rna', 'library_pcr_forward', 'library_pcr_reverse', 'library_pcr_polymerase',
                     'library_pcr', 'library_page_selection', 'library_page', 'library_library'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _runon_cap_graph(DEFINITION, FINAL_LIBRARY)
