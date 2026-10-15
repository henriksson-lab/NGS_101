from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, RLoopCapture, _part, _role_library, _role_library_scene
DEFINITION = Definition("DRIPc-seq — strand-specific sequencing of R-loop RNA",
        "10.1016/j.molcel.2016.05.032", "Sanz et al., Molecular Cell (2016)",
        "After S9.6 capture, DNase releases the R-loop RNA for reverse transcription; dUTP second-strand synthesis and UNG make the final library strand-specific.",
        ("R-loop RNA-derived cDNA",), "Single-end reads derive from the retained cDNA strand and therefore identify the RNA strand in the original hybrid.", "01_dripc_seq.html",
        caveat="The paper and detailed protocol specify indexed Illumina construction but do not print complete adapter sequences; adapter regions are shown by role.")

FINAL_LIBRARY, SEQ_PRIMERS = _role_library(
    DEFINITION.title + " final library", DEFINITION.payload), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Detailed chemistry was checked against the authors’ later protocol.")
SPEC.CAVEAT = ("The authors’ detailed protocol prints the P5/P7 library PCR primers and "
               "example adapter-trimming motifs. It names TruSeq Single index ligation "
               "adapters but does not print the selected index oligos or their end geometry.")
SPEC.oligos = lambda: (
    "<b>Library PCR P5</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′",
    "<b>Library PCR P7</b> — 5′-CAAGCAGAAGACGGCATACGAGAT-3′",
)

def _hybrid_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed graph for the DNA:RNA-hybrid capture protocols in this module."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cap = RLoopCapture('S9.6 antibody', 'RNA')
    chromatin = graph.input('rloop_chromatin', 'chromatin containing R-loops', MaterialKind.BIOLOGICAL, state=cap.rloop(), primary=True)
    extraction = graph.input('gentle_extraction', 'gentle nucleic-acid extraction', MaterialKind.CONDITION)
    extracted = graph.product('extracted_nucleic_acid', 'gently extracted R-loop-containing nucleic acid', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('extract', 'Gently extract nucleic acid', StepKind.EXTRACTION, ((chromatin, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (extracted,))
    restriction = graph.input('restriction_cocktail', 'R-loop-preserving restriction-enzyme cocktail', MaterialKind.CATALYST)
    digested = graph.product('restriction_fragments', 'restriction-fragmented R-loop nucleic acid', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('restriction_digest', 'Restriction-digest genomic DNA', StepKind.DIGESTION, ((extracted, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    antibody = graph.input('s96_antibody', 'S9.6 antibody', MaterialKind.SELECTION)
    beads = graph.input('capture_beads', 'antibody-capture beads', MaterialKind.SELECTION)
    captured = graph.product('captured_rloops', 'S9.6-captured R-loops', MaterialKind.MOLECULAR, state=cap.captured_component())
    rejected = graph.product('uncaptured_nucleic_acid', 'uncaptured nucleic acid', MaterialKind.WASTE)
    graph.react('immunoprecipitate', 'Immunoprecipitate RNA:DNA hybrids', StepKind.AFFINITY_CAPTURE, ((digested, InputRole.TEMPLATE), (antibody, InputRole.SELECTION), (beads, InputRole.SELECTION)), (captured, rejected))
    graph.discard(rejected, 'not retained by S9.6 immunoprecipitation')
    reviewed = ['rloop_chromatin', 'gentle_extraction', 'extract', 'restriction_cocktail', 'restriction_digest', 's96_antibody', 'capture_beads', 'immunoprecipitate']
    dnase = graph.input('dnase', 'DNase', MaterialKind.CATALYST)
    retained_rna = graph.product('retained_rloop_rna', 'retained R-loop RNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    dna_waste = graph.product('digested_dna', 'digested DNA component', MaterialKind.WASTE)
    graph.react('digest_dna', 'Digest the captured DNA component', StepKind.DIGESTION, ((captured, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (retained_rna, dna_waste))
    graph.discard(dna_waste, 'DNase removes the DNA component before RNA copying')
    rt_primers = graph.input('rt_primer_mix', 'poly(dT) and random-hexamer RT primers', MaterialKind.OLIGO)
    rtase = graph.input('reverse_transcriptase', 'reverse transcriptase', MaterialKind.CATALYST)
    cdna = graph.product('first_strand_cdna', 'R-loop RNA-derived cDNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    graph.react('reverse_transcription', 'Reverse-transcribe retained R-loop RNA', StepKind.REVERSE_TRANSCRIPTION, ((retained_rna, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (cdna,))
    second_primer = graph.input('second_strand_primer', 'second-strand primer', MaterialKind.OLIGO)
    second_polymerase = graph.input('second_strand_polymerase', 'dUTP second-strand polymerase', MaterialKind.CATALYST)
    dutp = graph.input('dutp', 'dUTP-containing nucleotide mix', MaterialKind.CONDITION)
    duplex = graph.product('dutp_cdna_duplex', 'dUTP-marked cDNA duplex', MaterialKind.MOLECULAR, state=_role_library_scene('dUTP-marked cDNA', ('R-loop RNA-derived cDNA', 'dUTP strand'), ligated=False))
    graph.react('second_strand', 'Synthesize the dUTP-marked second strand', StepKind.EXTENSION, ((cdna, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL), (second_polymerase, InputRole.CATALYST), (dutp, InputRole.CONDITION)), (duplex,))
    sonication = graph.input('sonication', 'cDNA sonication condition', MaterialKind.CONDITION)
    fragmented = graph.product('fragmented_cdna', 'fragmented dUTP-marked cDNA', MaterialKind.MOLECULAR, state=duplex.state)
    graph.react('fragment_cdna', 'Sonicate the cDNA duplex', StepKind.FRAGMENTATION, ((duplex, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragmented,))
    adapters = graph.input('indexed_adapters', 'indexed library adapters', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    tagged = graph.product('adapter_ligated_cdna', 'adapter-ligated cDNA', MaterialKind.MOLECULAR, state=fragmented.state, inferred=True)
    graph.react('adapter_ligation', 'Ligate indexed adapters', StepKind.LIGATION, ((fragmented, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=AdapterLigationGeometry.unavailable('The source specifies indexed adapters but does not print their ligation-end geometry.'))
    ung = graph.input('ung', 'uracil-DNA glycosylase', MaterialKind.CATALYST)
    retained = graph.product('retained_cdna_strand', 'retained cDNA library strand', MaterialKind.MOLECULAR, state=tagged.state, inferred=True)
    uracil_waste = graph.product('removed_dutp_strand', 'removed dUTP-marked strand', MaterialKind.WASTE)
    graph.react('remove_dutp_strand', 'Remove the dUTP-containing strand', StepKind.DIGESTION, ((tagged, InputRole.TEMPLATE), (ung, InputRole.CATALYST)), (retained, uracil_waste))
    graph.discard(uracil_waste, 'UNG removes the dUTP-containing strand before PCR')
    pcr_input = retained
    reviewed.extend(('dnase', 'digest_dna', 'rt_primer_mix', 'reverse_transcriptase', 'reverse_transcription', 'second_strand_primer', 'second_strand_polymerase', 'dutp', 'second_strand', 'sonication', 'fragment_cdna', 'indexed_adapters', 'dna_ligase', 'adapter_ligation', 'ung', 'remove_dutp_strand'))
    forward = graph.input('library_forward_primer', 'library PCR forward primer', MaterialKind.OLIGO, inferred=True)
    reverse = graph.input('library_reverse_primer', 'library PCR reverse primer', MaterialKind.OLIGO, inferred=True)
    pcr_polymerase = graph.input('pcr_polymerase', 'library PCR polymerase', MaterialKind.CATALYST)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react('library_pcr', 'PCR-amplify the sequencing library', StepKind.PCR, ((pcr_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(pcr_input.state if isinstance(pcr_input.state, Construct) else final, reason='The authoritative source does not expose one exact selected library-primer pair.'))
    graph.sequence(library, SequencingPlan(limitation=definition.caveat or definition.sequencing))
    reviewed.extend(('library_forward_primer', 'library_reverse_primer', 'pcr_polymerase', 'library_pcr', 'sequencing_library'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _hybrid_graph(d, final)
SPEC.protocol_graph = protocol_graph
