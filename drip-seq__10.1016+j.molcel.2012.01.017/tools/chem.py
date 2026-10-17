from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, RLoopCapture, _part, _role_library, _role_library_scene
DEFINITION = Definition("DRIP-seq — DNA:RNA immunoprecipitation sequencing",
        "10.1016/j.molcel.2012.01.017", "Ginno et al., Molecular Cell (2012)",
        "Gentle extraction and restriction digestion preserve R-loops for S9.6 immunoprecipitation; a conventional dsDNA library sequences the captured DNA component.",
        ("captured R-loop DNA",), "Illumina GAIIx and HiSeq reads report captured restriction-fragment DNA without strand identity; the defining Methods do not state read mode or length.", "01_drip_seq.html",
        caveat="The 2012 Methods state only that the S9.6 pull-down, the RNase H-pre-treated control and input DNA were sonicated, size-selected and ligated to Illumina barcoded adaptors. End repair, dA-tailing and library PCR are drawn from the authors' later (2019) DRIP-seq library branch and marked inferred; adapter and primer sequences are shown by role.")
LATER_PROTOCOL = "Sanz & Chédin, Nat Protoc (2019), doi:10.1038/s41596-019-0159-1, Procedure steps 44–51 (DRIP-seq skips step 50)"

FINAL_LIBRARY, SEQ_PRIMERS = _role_library(
    DEFINITION.title + " final library", DEFINITION.payload), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Library steps not stated in 2012 follow the authors’ later protocol (Sanz &amp; Chédin, <i>Nat Protoc</i> 2019) and are marked inferred.")

def _hybrid_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed graph for the DNA:RNA-hybrid capture protocols in this module."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cap = RLoopCapture('S9.6 antibody', 'DNA')
    chromatin = graph.input('rloop_chromatin', 'chromatin containing R-loops', MaterialKind.BIOLOGICAL, state=cap.rloop(), primary=True)
    extraction = graph.input('gentle_extraction', 'SDS/proteinase K at 37 °C, phenol-chloroform, ethanol precipitation', MaterialKind.CONDITION)
    extracted = graph.product('extracted_nucleic_acid', 'gently extracted R-loop-containing nucleic acid', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('extract', 'Gently extract nucleic acid', StepKind.EXTRACTION, ((chromatin, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (extracted,))
    restriction = graph.input('restriction_cocktail', 'HindIII, EcoRI, BsrGI, XbaI and SspI', MaterialKind.CATALYST)
    digested = graph.product('restriction_fragments', 'restriction-fragmented R-loop nucleic acid', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('restriction_digest', 'Restriction-digest genomic DNA', StepKind.DIGESTION, ((extracted, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    antibody = graph.input('s96_antibody', 'S9.6 antibody', MaterialKind.SELECTION)
    beads = graph.input('capture_beads', 'antibody-capture beads', MaterialKind.SELECTION)
    captured = graph.product('captured_rloops', 'S9.6-captured R-loops', MaterialKind.MOLECULAR, state=cap.captured_component())
    rejected = graph.product('uncaptured_nucleic_acid', 'uncaptured nucleic acid', MaterialKind.WASTE)
    graph.react('immunoprecipitate', 'Immunoprecipitate RNA:DNA hybrids (MeDIP-style, S9.6, no denaturation)', StepKind.AFFINITY_CAPTURE, ((digested, InputRole.TEMPLATE), (antibody, InputRole.SELECTION), (beads, InputRole.SELECTION)), (captured, rejected))
    graph.discard(rejected, 'not retained by S9.6 immunoprecipitation')
    reviewed = ['rloop_chromatin', 'gentle_extraction', 'extract', 'restriction_cocktail', 'restriction_digest', 's96_antibody', 'capture_beads', 'immunoprecipitate']
    sonication = graph.input('sonication', 'post-capture sonication condition', MaterialKind.CONDITION)
    fragmented = graph.product('captured_dna_fragments', 'sonicated captured DNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    graph.react('fragment_captured_dna', 'Sonicate captured DNA', StepKind.FRAGMENTATION, ((captured, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragmented,))
    size_window = graph.input('size_selection', 'size selection (window not stated)', MaterialKind.SELECTION)
    sized = graph.product('size_selected_dna', 'size-selected captured DNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    out_of_range = graph.product('out_of_range_fragments', 'fragments outside the selected size window', MaterialKind.WASTE)
    graph.react('select_size', 'Size-select sonicated DNA', StepKind.SELECTION, ((fragmented, InputRole.TEMPLATE), (size_window, InputRole.SELECTION)), (sized, out_of_range), changes_molecule=False)
    graph.discard(out_of_range, 'outside the size-selected window (2012 Methods; window not stated)')
    repair = graph.input('end_repair_enzymes', 'end-repair enzyme mix (NEBNext End Repair in the 2019 protocol)', MaterialKind.CATALYST, inferred=True, source=LATER_PROTOCOL)
    repaired = graph.product('end_repaired_dna', 'end-repaired captured DNA', MaterialKind.MOLECULAR, state=sized.state, inferred=True, source=LATER_PROTOCOL)
    graph.react('end_repair', 'INFERRED — End-repair captured DNA', StepKind.END_REPAIR, ((sized, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input('a_tailing_enzyme', 'dA-tailing enzyme', MaterialKind.CATALYST, inferred=True, source=LATER_PROTOCOL)
    tailed = graph.product('da_tailed_dna', 'dA-tailed captured DNA', MaterialKind.MOLECULAR, state=_role_library_scene('dA-tailed captured DNA', ('captured DNA fragment',), ligated=False), inferred=True, source=LATER_PROTOCOL)
    graph.react('da_tail', 'INFERRED — Add one 3-prime dA to each end', StepKind.TAILING, ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input('library_adapters', 'Illumina barcoded adaptors (sequences not printed)', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    tagged = graph.product('adapter_ligated_dna', 'adapter-ligated captured DNA', MaterialKind.MOLECULAR, state=_role_library_scene('non-strand-specific DRIP-seq library', ('captured DNA fragment',)), inferred=True)
    graph.react('adapter_ligation', 'Ligate adapters', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=AdapterLigationGeometry.unavailable('The 2012 Methods name Illumina barcoded adaptors without their sequences or end geometry; the 2019 protocol uses TruSeq single-index adapters after dA-tailing.'))
    pcr_input = tagged
    reviewed.extend(('sonication', 'fragment_captured_dna', 'size_selection', 'select_size', 'size_selected_dna', 'out_of_range_fragments', 'end_repair_enzymes', 'end_repair', 'a_tailing_enzyme', 'da_tail', 'library_adapters', 'dna_ligase', 'adapter_ligation'))
    forward = graph.input('library_forward_primer', 'library PCR forward primer', MaterialKind.OLIGO, inferred=True, source=LATER_PROTOCOL)
    reverse = graph.input('library_reverse_primer', 'library PCR reverse primer', MaterialKind.OLIGO, inferred=True, source=LATER_PROTOCOL)
    pcr_polymerase = graph.input('pcr_polymerase', 'library PCR polymerase', MaterialKind.CATALYST)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react('library_pcr', 'INFERRED — PCR-amplify the sequencing library', StepKind.PCR, ((pcr_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(pcr_input.state if isinstance(pcr_input.state, Construct) else final, reason='The authoritative source does not expose one exact selected library-primer pair.'))
    graph.sequence(library, SequencingPlan(limitation=definition.caveat or definition.sequencing))
    reviewed.extend(('library_forward_primer', 'library_reverse_primer', 'pcr_polymerase', 'library_pcr', 'sequencing_library'))
    graph.annotate(Highlight('immunoprecipitate', 'Controls are separate libraries',
                             'The same digest pre-treated overnight with RNase H is immunoprecipitated in parallel, and input DNA is also sequenced; DRIP-seq peaks require >10-fold enrichment over both (2012 Methods).'))
    graph.annotate(Highlight('adapter_ligation', '2012 versus 2019 library construction',
                             'The 2012 Methods name only sonication, size selection and ligation to Illumina barcoded adaptors. End repair, dA-tailing, TruSeq single-index adapters and P5/P7 PCR belong to the authors’ 2019 protocol and are not assigned to the 2012 libraries.'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _hybrid_graph(d, final)
SPEC.protocol_graph = protocol_graph
