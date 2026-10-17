from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, RLoopCapture, _part, _role_library, _role_library_scene
DEFINITION = Definition("DRIPc-seq — strand-specific sequencing of R-loop RNA",
        "10.1016/j.molcel.2016.05.032", "Sanz et al., Molecular Cell (2016)",
        "After S9.6 capture, DNase releases the R-loop RNA for reverse transcription; dUTP second-strand synthesis and UNG make the final library strand-specific.",
        ("R-loop RNA-derived cDNA",), "Reads derive from the retained cDNA strand and therefore identify the RNA strand in the original hybrid. The 2016 article reports Illumina HiSeq 2000 sequencing without a read length; the 2019 protocol mentions single-end 100 bp only as an example case.", "01_dripc_seq.html",
        caveat="The paper and detailed protocol specify indexed Illumina construction but do not print complete adapter sequences; adapter regions are shown by role.")

FINAL_LIBRARY, SEQ_PRIMERS = _role_library(
    DEFINITION.title + " final library", DEFINITION.payload), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Library steps (end repair, dA-tailing, ligation, UNG, PCR, size selection) follow the authors’ later protocol, <a href=\"https://doi.org/10.1038/s41596-019-0159-1\">Sanz &amp; Chédin, Nat Protoc (2019)</a>; the 2016 article states only sonication, library construction, a UDG step before PCR and HiSeq 2000 sequencing.")
SPEC.CAVEAT = ("The authors’ detailed protocol prints the P5/P7 library PCR primers and "
               "example adapter-trimming motifs. It names TruSeq Single index ligation "
               "adapters but does not print the selected index oligos or their end geometry.")
# Sanz & Chédin 2019 Nat Protoc Reagents (PMC6615061 text lines617-621), verbatim.
PCR_PRIMER_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA"
PCR_PRIMER_P7 = "CAAGCAGAAGACGGCATACGAGAT"


def pcr_primer_p5() -> Construct:
    return Construct([Segment("PCR primer 1.0 P5", PCR_PRIMER_P5, "p5")], name="PCR primer 1.0 P5")


def pcr_primer_p7() -> Construct:
    return Construct([Segment("PCR primer 2.0 P7", PCR_PRIMER_P7, "p7")], name="PCR primer 2.0 P7")


SPEC.oligos = lambda: (
    f"<b>PCR primer 1.0 P5</b> — 5′-{PCR_PRIMER_P5}-3′",
    f"<b>PCR primer 2.0 P7</b> — 5′-{PCR_PRIMER_P7}-3′",
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
    dnase = graph.input('dnase', 'DNase I (NEB); EDTA/75 °C inactivation in the 2019 protocol', MaterialKind.CATALYST)
    retained_rna = graph.product('retained_rloop_rna', 'retained R-loop RNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    dna_waste = graph.product('digested_dna', 'digested DNA component', MaterialKind.WASTE)
    graph.react('digest_dna', 'Digest the captured DNA component', StepKind.DIGESTION, ((captured, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (retained_rna, dna_waste))
    graph.discard(dna_waste, 'DNase I degrades the DNA component; RNA is then ethanol-precipitated')
    rt_primers = graph.input('rt_primer_mix', 'poly(T) and random-hexamer primers (in iScript RT supermix)', MaterialKind.OLIGO)
    rtase = graph.input('reverse_transcriptase', 'iScript reverse transcriptase (Bio-Rad)', MaterialKind.CATALYST)
    cdna = graph.product('first_strand_cdna', 'R-loop RNA-derived cDNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    graph.react('reverse_transcription', 'Reverse-transcribe retained R-loop RNA', StepKind.REVERSE_TRANSCRIPTION, ((retained_rna, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (cdna,))
    # 2019 protocol steps 39-40: RNase H nicking with E. coli DNA ligase and DNA polymerase I
    # (no added second-strand oligo); the 2016 article states only dUTP replacing dTTP.
    second_rnaseh = graph.input('second_strand_rnase_h', 'RNase H (nicks the RNA in the RNA:cDNA hybrid)', MaterialKind.CATALYST)
    second_ligase = graph.input('second_strand_ligase', 'E. coli DNA ligase (NAD+)', MaterialKind.CATALYST)
    second_polymerase = graph.input('second_strand_polymerase', 'E. coli DNA polymerase I', MaterialKind.CATALYST)
    dutp = graph.input('dutp', 'dUTP-containing nucleotide mix', MaterialKind.CONDITION)
    duplex = graph.product('dutp_cdna_duplex', 'dUTP-marked cDNA duplex', MaterialKind.MOLECULAR, state=_role_library_scene('dUTP-marked cDNA', ('R-loop RNA-derived cDNA', 'dUTP strand'), ligated=False))
    graph.react('second_strand', 'Synthesize the dUTP-marked second strand', StepKind.EXTENSION, ((cdna, InputRole.TEMPLATE), (second_rnaseh, InputRole.CATALYST), (second_ligase, InputRole.CATALYST), (second_polymerase, InputRole.CATALYST), (dutp, InputRole.CONDITION)), (duplex,))
    sonication = graph.input('sonication', 'Bioruptor sonication to ~200 bp (2016) / 200–300 bp average (2019)', MaterialKind.CONDITION)
    fragmented = graph.product('fragmented_cdna', 'fragmented dUTP-marked cDNA', MaterialKind.MOLECULAR, state=duplex.state)
    graph.react('fragment_cdna', 'Sonicate the cDNA duplex', StepKind.FRAGMENTATION, ((duplex, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragmented,))
    end_repair = graph.input('end_repair_module', 'NEBNext End Repair module', MaterialKind.CATALYST)
    repaired = graph.product('end_repaired_cdna', 'end-repaired cDNA fragments', MaterialKind.MOLECULAR, state=fragmented.state)
    graph.react('end_repair', 'End-repair sonicated cDNA (2019 step 44)', StepKind.END_REPAIR, ((fragmented, InputRole.TEMPLATE), (end_repair, InputRole.CATALYST)), (repaired,))
    klenow = graph.input('klenow_exo', 'Klenow fragment (3′→5′ exo−) with dATP', MaterialKind.CATALYST)
    a_tailed = graph.product('a_tailed_cdna', 'dA-tailed cDNA fragments', MaterialKind.MOLECULAR, state=fragmented.state)
    graph.react('a_tailing', 'dA-tail the repaired cDNA (2019 step 46)', StepKind.TAILING, ((repaired, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (a_tailed,))
    adapters = graph.input('indexed_adapters', 'Illumina TruSeq Single index adapters (selected index not printed)', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('dna_ligase', 'NEB Quick ligase', MaterialKind.CATALYST)
    tagged = graph.product('adapter_ligated_cdna', 'adapter-ligated cDNA', MaterialKind.MOLECULAR, state=fragmented.state, inferred=True)
    graph.react('adapter_ligation', 'Ligate indexed adapters', StepKind.LIGATION, ((a_tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=AdapterLigationGeometry.unavailable('The source specifies indexed adapters but does not print their ligation-end geometry.'))
    ung = graph.input('ung', 'AmpErase uracil N-glycosylase', MaterialKind.CATALYST)
    retained = graph.product('retained_cdna_strand', 'retained cDNA library strand', MaterialKind.MOLECULAR, state=tagged.state, inferred=True)
    uracil_waste = graph.product('removed_dutp_strand', 'removed dUTP-marked strand', MaterialKind.WASTE)
    graph.react('remove_dutp_strand', 'Remove the dUTP-containing strand', StepKind.DIGESTION, ((tagged, InputRole.TEMPLATE), (ung, InputRole.CATALYST)), (retained, uracil_waste))
    graph.discard(uracil_waste, 'UNG removes the dUTP-containing strand before PCR')
    pcr_input = retained
    reviewed.extend(('dnase', 'digest_dna', 'rt_primer_mix', 'reverse_transcriptase', 'reverse_transcription', 'second_strand_rnase_h', 'second_strand_ligase', 'second_strand_polymerase', 'dutp', 'second_strand', 'sonication', 'fragment_cdna', 'end_repair_module', 'end_repair', 'klenow_exo', 'a_tailing', 'indexed_adapters', 'dna_ligase', 'adapter_ligation', 'ung', 'remove_dutp_strand'))
    forward = graph.input('library_forward_primer', 'PCR primer 1.0 P5 (sequence printed in the 2019 protocol)', MaterialKind.OLIGO, state=pcr_primer_p5())
    reverse = graph.input('library_reverse_primer', 'PCR primer 2.0 P7 (sequence printed in the 2019 protocol)', MaterialKind.OLIGO, state=pcr_primer_p7())
    pcr_polymerase = graph.input('pcr_polymerase', 'Phusion Flash master mix (10–15 cycles)', MaterialKind.CATALYST)
    amplified = graph.product('amplified_library', 'amplified library', MaterialKind.MOLECULAR, state=final, inferred=True)
    graph.react('library_pcr', 'PCR-amplify the sequencing library', StepKind.PCR, ((pcr_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(pcr_input.state if isinstance(pcr_input.state, Construct) else final, reason='The authoritative source does not expose one exact selected library-primer pair.'))
    beads = graph.input('ampure_size_selection', 'two-ratio AMPure XP selection (0.65× then added beads)', MaterialKind.SELECTION)
    library = graph.product('sequencing_library', definition.title + ' sequencing library (200–500 bp)', MaterialKind.LIBRARY, state=final, inferred=True)
    large = graph.product('oversized_fragments', 'fragments above ~500 bp', MaterialKind.WASTE)
    small = graph.product('undersized_fragments', 'fragments below ~200 bp in the supernatant', MaterialKind.WASTE)
    graph.react('size_select', 'Two-step AMPure size selection (2019 steps 53–56)', StepKind.SELECTION, ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (library, large, small))
    graph.discard(large, 'bound to the first 0.65× AMPure beads, which are discarded')
    graph.discard(small, 'left in the second-step supernatant')
    graph.annotate(SourceDiscrepancy('digest_dna', 'DNase I conditions differ between the article and the later protocol',
        (SourceClaim('Sanz et al. 2016 Experimental Procedures, PMC4955522 text line137', '4 U DNase I, 35 min at 37 °C'),
         SourceClaim('Sanz & Chédin 2019 Nat Protoc step 26, PMC6615061 text lines1020–1036', '6 U DNase I (3 µL) in 100 µL, 45 min at 37 °C, then EDTA and 75 °C 15 min')),
        'Leave the DNase dose unassigned in the graph; the operation (complete DNA digestion) is shared.',
        'Both values are printed in the respective originals.',
        'The 2019 protocol is a later, detailed procedure from the same laboratory; whether the change is a deliberate revision is not stated.'))
    graph.sequence(library, SequencingPlan(limitation=definition.caveat or definition.sequencing))
    reviewed.extend(('library_forward_primer', 'library_reverse_primer', 'pcr_polymerase', 'library_pcr', 'amplified_library', 'ampure_size_selection', 'size_select', 'sequencing_library', 'oversized_fragments', 'undersized_fragments'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _hybrid_graph(d, final)
SPEC.protocol_graph = protocol_graph
