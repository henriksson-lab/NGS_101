from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import illumina as il
from chemdraw import Segment
from endprep import dA_tailed_scene, repair_and_dA_tail
from hybrid_base_protocols import Definition, RLoopCapture, _part, _role_library_scene, make_spec
def rchip_adapter_geometry() -> AdapterLigationGeometry:
    """Authors' 2019 R-ChIP protocol, Table 1: one-ended adapter geometry."""
    sc = Scene()
    sc.strand("oligo A", [Segment("A paired core", "GATCGGAAGAGC")],
              label="5′-phosphorylated adapter oligo A", mod5="p")
    sc.anneal("oligo B", [Segment("B paired core", "GCTCTTCCGATC"),
                            Segment("B 3-prime T overhang", "T")],
              to="oligo A", pair=("B paired core", "A paired core"),
              label="adapter oligo B")
    sc.note("oligo A", "unpaired 3′ tail: GTCGTGTAGGGAAAGAGTGT")
    sc.note("oligo B", "unpaired 5′ tail: AGACGTGT")
    sc.ligation_end("oligo B", "B 3-prime T overhang", kind="cohesive",
                    polarity="3-prime", owner="adapter", text="3′-T ligating end")
    substrate = dA_tailed_scene(repair_and_dA_tail(Construct([
        Segment("R-loop DNA-derived insert", "X" * 24, placeholder=True)])),
        label="unobstructed dA-tailed DNA end")
    return AdapterLigationGeometry.known(
        "cohesive", substrate, sc,
        statement="Table 1 of the authors’ 2019 protocol prints both adapter strands; "
                  "their 12-base antiparallel overlap leaves a 3′-T end compatible "
                  "with the dA-tailed captured DNA.")


def _rchip_library() -> Construct:
    return Construct([
        Segment("tailed N9 handle", "CAAGCAGAAGACGGCATACGAG", "p7"),
        _part("N9 random-primer remnant", "N", 9),
        _part("R-loop DNA-derived insert"),
        _part("one ligated Illumina adapter", "R", 18),
    ], name="strand-specific R-ChIP library")


def _rchip_library_scene() -> Scene:
    sc = Scene.duplex(list(_rchip_library()), label="strand-specific R-ChIP library")
    sc.junction("top", "R-loop DNA-derived insert", "one ligated Illumina adapter",
                "adapter ligation")
    return sc


DEFINITION = Definition("R-ChIP — in-vivo R-loop capture with inactive RNASEH1",
        "10.1016/j.molcel.2017.10.008", "Chen et al., Molecular Cell (2017)",
        "V5-tagged RNASEH1 D210N binds R-loops in fixed cells; a tailed N9 primer and one-ended adapter ligation produce a strand-specific DNA-component library.",
        ("tailed N9 end", "R-loop DNA-derived insert"), "Single-end Illumina reads preserve the orientation of the DNA strand captured from each R-loop.", "01_r_chip.html",
        caveat="The authors’ detailed 2019 protocol prints both one-ended adapter strands and both PCR primers; the final amplified adapter shell remains partly inferred.")

FINAL_LIBRARY, SEQ_PRIMERS = _rchip_library(), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Detailed chemistry was checked against the authors’ later protocol.")
SPEC.CAVEAT = ("The defining paper prints the tailed N9 primer; the authors’ detailed "
               "2019 protocol also prints both adapter strands and the PCR primers. "
               "The final adapter body remains a role-level segment pending reconstruction "
               "of its one-ended amplification product.")
SPEC.SEQUENCING_UNAVAILABLE = (
    "The library-preparation PCR primers are published; a separate instrument "
    "sequencing-primer synthesis oligo is not specified for this run.")
SPEC.oligos = lambda: (
    "<b>Tailed N9 primer</b> — 5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′",
    "<b>Adapter oligo A</b> — 5′-/Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3′",
    "<b>Adapter oligo B</b> — 5′-AGACGTGTGCTCTTCCGATCT-3′",
    "<b>PCR primer</b> — 5′-CAAGCAGAAGACGGCATACGAG-3′",
    "<b>Barcoded PCR primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
)


def adapter_geometry():
    return rchip_adapter_geometry()


def _hybrid_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed graph for the DNA:RNA-hybrid capture protocols in this module."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cap = RLoopCapture('V5-tagged catalytically inactive RNASEH1 D210N', 'DNA', in_vivo=True)
    chromatin = graph.input('rloop_chromatin', 'chromatin containing R-loops', MaterialKind.BIOLOGICAL, state=cap.rloop(), primary=True)
    formaldehyde = graph.input('formaldehyde', 'formaldehyde cross-linking reagent', MaterialKind.CONDITION)
    fixed = graph.product('fixed_chromatin', 'cross-linked R-loop chromatin', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('crosslink', 'Cross-link cells', StepKind.CROSS_LINKING, ((chromatin, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    sonication = graph.input('sonication', 'chromatin sonication condition', MaterialKind.CONDITION)
    sheared = graph.product('sheared_chromatin', 'sheared fixed chromatin', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('shear', 'Shear cross-linked chromatin', StepKind.FRAGMENTATION, ((fixed, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (sheared,))
    v5_antibody = graph.input('v5_antibody', 'anti-V5 antibody', MaterialKind.SELECTION)
    inactive_rnaseh = graph.input('inactive_rnaseh1', 'V5–RNASEH1 D210N R-loop binder', MaterialKind.SELECTION)
    captured = graph.product('captured_rloop_dna', 'captured R-loop chromatin', MaterialKind.MOLECULAR, state=cap.captured_component())
    rejected = graph.product('uncaptured_chromatin', 'uncaptured chromatin', MaterialKind.WASTE)
    graph.react('immunoprecipitate', 'Immunoprecipitate V5–RNASEH1-bound R-loops', StepKind.AFFINITY_CAPTURE, ((sheared, InputRole.TEMPLATE), (inactive_rnaseh, InputRole.SELECTION), (v5_antibody, InputRole.SELECTION)), (captured, rejected))
    graph.discard(rejected, 'not retained by anti-V5 immunoprecipitation')
    reviewed = ['rloop_chromatin', 'formaldehyde', 'crosslink', 'sonication', 'shear', 'v5_antibody', 'inactive_rnaseh1', 'immunoprecipitate']
    denaturing = graph.input('denaturing_condition', 'captured-DNA denaturation', MaterialKind.CONDITION)
    single = graph.product('single_stranded_captured_dna', 'single-stranded captured R-loop DNA', MaterialKind.MOLECULAR, state=cap.captured_component())
    partner = graph.product('released_complement', 'released complementary strand', MaterialKind.WASTE)
    graph.react('denature', 'Denature captured R-loop DNA', StepKind.STRAND_SEPARATION, ((captured, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)), (single, partner))
    graph.discard(partner, 'the copied captured strand defines the directional library')
    n9 = graph.input('tailed_n9_primer', 'tailed N9 random primer', MaterialKind.OLIGO)
    polymerase = graph.input('copy_polymerase', 'DNA polymerase', MaterialKind.CATALYST)
    primed_scene = Scene()
    primed_scene.strand('copied molecule', [_part('N9-primer handle'), _part('N9 random primer', 'N', 9), _part('R-loop DNA-derived copy')], label='unidirectional R-ChIP copy')
    primed = graph.product('n9_copied_dna', 'tailed-N9 copied R-loop DNA', MaterialKind.MOLECULAR, state=primed_scene)
    graph.react('random_primer_copy', 'Copy captured DNA with the tailed N9 primer', StepKind.EXTENSION, ((single, InputRole.TEMPLATE), (n9, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (primed,))
    tailing = graph.input('a_tailing_enzyme', 'dA-tailing enzyme', MaterialKind.CATALYST)
    tailed = graph.product('one_end_da_tailed', 'one-end dA-tailed copied DNA', MaterialKind.MOLECULAR, state=primed_scene)
    graph.react('da_tail', 'Add one 3-prime dA to the unobstructed end', StepKind.TAILING, ((primed, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input('illumina_adapter', 'one source-written Illumina adapter', MaterialKind.ADAPTER)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    tagged = graph.product('one_end_ligated_dna', 'one-end-adapter-ligated DNA', MaterialKind.MOLECULAR, state=_rchip_library_scene(), inferred=True)
    graph.react('adapter_ligation', 'Ligate one Illumina adapter', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=rchip_adapter_geometry())
    pcr_input = tagged
    reviewed.extend(('denaturing_condition', 'denature', 'tailed_n9_primer', 'copy_polymerase', 'random_primer_copy', 'a_tailing_enzyme', 'da_tail', 'illumina_adapter', 'dna_ligase', 'adapter_ligation'))
    forward = graph.input('library_forward_primer', 'library PCR forward primer', MaterialKind.OLIGO, inferred=False)
    reverse = graph.input('library_reverse_primer', 'library PCR reverse primer', MaterialKind.OLIGO, inferred=False)
    pcr_polymerase = graph.input('pcr_polymerase', 'library PCR polymerase', MaterialKind.CATALYST)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react('library_pcr', 'PCR-amplify the sequencing library', StepKind.PCR, ((pcr_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(pcr_input.state if isinstance(pcr_input.state, Construct) else final, reason='The published R-ChIP PCR primers cannot be placed on the still inferred final amplification shell.'))
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
