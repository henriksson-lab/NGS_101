from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, RLoopCapture, _part, _role_library, _role_library_scene
DEFINITION = Definition("DRIP-seq — DNA:RNA immunoprecipitation sequencing",
        "10.1016/j.molcel.2012.01.017", "Ginno et al., Molecular Cell (2012)",
        "Gentle extraction and restriction digestion preserve R-loops for S9.6 immunoprecipitation; a conventional dsDNA library sequences the captured DNA component.",
        ("captured R-loop DNA",), "Single-end Illumina reads report captured restriction fragments without strand identity.", "01_drip_seq.html",
        caveat="The defining paper names the library reagents but does not expose complete adapter or sequencing-primer sequences; adapter regions are shown by role.")

FINAL_LIBRARY, SEQ_PRIMERS = _role_library(
    DEFINITION.title + " final library", DEFINITION.payload), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Detailed chemistry was checked against the authors’ later protocol.")

def _hybrid_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed graph for the DNA:RNA-hybrid capture protocols in this module."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cap = RLoopCapture('S9.6 antibody', 'DNA')
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
    sonication = graph.input('sonication', 'post-capture sonication condition', MaterialKind.CONDITION)
    fragmented = graph.product('captured_dna_fragments', 'sonicated captured DNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    graph.react('fragment_captured_dna', 'Sonicate captured DNA', StepKind.FRAGMENTATION, ((captured, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragmented,))
    repair = graph.input('end_repair_enzymes', 'end-repair enzyme mix', MaterialKind.CATALYST)
    repaired = graph.product('end_repaired_dna', 'end-repaired captured DNA', MaterialKind.MOLECULAR, state=fragmented.state)
    graph.react('end_repair', 'End-repair captured DNA', StepKind.END_REPAIR, ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input('a_tailing_enzyme', 'dA-tailing enzyme', MaterialKind.CATALYST)
    tailed = graph.product('da_tailed_dna', 'dA-tailed captured DNA', MaterialKind.MOLECULAR, state=_role_library_scene('dA-tailed captured DNA', ('captured DNA fragment',), ligated=False))
    graph.react('da_tail', 'Add one 3-prime dA to each end', StepKind.TAILING, ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input('library_adapters', 'standard dsDNA library adapters', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    tagged = graph.product('adapter_ligated_dna', 'adapter-ligated captured DNA', MaterialKind.MOLECULAR, state=_role_library_scene('non-strand-specific DRIP-seq library', ('captured DNA fragment',)), inferred=True)
    graph.react('adapter_ligation', 'Ligate adapters', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=AdapterLigationGeometry.unavailable('The detailed source specifies a standard dsDNA library but does not print the adapter-end oligos.'))
    pcr_input = tagged
    reviewed.extend(('sonication', 'fragment_captured_dna', 'end_repair_enzymes', 'end_repair', 'a_tailing_enzyme', 'da_tail', 'library_adapters', 'dna_ligase', 'adapter_ligation'))
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
