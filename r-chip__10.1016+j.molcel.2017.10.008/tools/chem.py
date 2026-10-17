from adapter_ligation import AdapterLigationGeometry
from dataclasses import replace
from chemdraw import Construct, Scene, complement_segments, feature, revcomp
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind
import illumina as il
from chemdraw import Segment
from hybrid_base_protocols import Definition, RLoopCapture, _part, make_spec
# Chen et al., Nat. Protoc. 2019 (PMC6604627), Table 1 — transcribed as printed.
N9_HANDLE = "CAAGCAGAAGACGGCATACGAG"          # 5' tail of the /invddt/ N9 primer
OLIGO_A = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # 5'-/Phos/
OLIGO_B = "AGACGTGTGCTCTTCCGATCT"
PCR_PRIMER = "CAAGCAGAAGACGGCATACGAG"
BARCODE_PRIMER = "AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"
N9_PRIMER = N9_HANDLE + "N" * 9
# The barcode primer is P5, a 5-nt flexible barcode and the TruSeq Read 1 sequence.
assert BARCODE_PRIMER == il.P5 + "NNNNN" + il.TRUSEQ_READ1
SAMPLE_BARCODE = feature("rchip_barcode_primer_nnnnn", "sample_index", "unknown")


def _a_core() -> Segment:
    return Segment("oligo A paired core", OLIGO_A[:12])


def _a_tail() -> Segment:
    return Segment("oligo A 3-prime tail", OLIGO_A[12:])


def _copy_segments() -> list[Segment]:
    return [Segment("tailed-N9 handle", N9_HANDLE, "p7"),
            _part("N9 random-primer remnant", "N", 9),
            _part("R-loop DNA-derived insert")]


def _adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("oligo A", [_a_core(), _a_tail()],
              label="5′-phosphorylated adapter oligo A", mod5="p")
    sc.anneal("oligo B", [Segment("B 5-prime tail", OLIGO_B[:8]),
                            Segment("B paired core", OLIGO_B[8:20]),
                            Segment("B 3-prime T overhang", OLIGO_B[20:])],
              to="oligo A", pair=("B paired core", "oligo A paired core"),
              label="adapter oligo B", unpaired=("B 5-prime tail",))
    sc.ligation_end("oligo B", "B 3-prime T overhang", kind="cohesive",
                    polarity="3-prime", owner="adapter", text="3′-T ligating end")
    return sc


def _one_end_tailed_scene() -> Scene:
    """Copied DNA: 5′ handle overhang at one end, a single 3′ dA at the other."""
    copy = _copy_segments()
    sc = Scene()
    sc.strand("copy", [*copy, Segment("3-prime dA", "A")],
              label="tailed-N9 copy strand", mod5="invddT")
    template = complement_segments(copy[1:])
    sc.anneal("template", template, to="copy",
              pair=("R-loop DNA-derived insert'", "R-loop DNA-derived insert"),
              label="captured template DNA")
    sc.overhang("copy", "3-prime dA", polarity="3-prime", owner="substrate", text="3′ dA")
    sc.note("copy", "5′ handle overhang blocks ligation at this end")
    return sc


def rchip_adapter_geometry() -> AdapterLigationGeometry:
    """Authors' 2019 R-ChIP protocol, Table 1 and steps 61–73: one-ended adapter geometry."""
    return AdapterLigationGeometry.known(
        "cohesive", _one_end_tailed_scene(), _adapter_scene(),
        statement="Table 1 of the authors’ 2019 protocol prints both adapter strands; "
                  "their 12-base antiparallel overlap leaves a 3′-T end compatible "
                  "with the single dA-tailed end of the copied DNA; the tailed-N9 "
                  "5′ overhang leaves the other end unligatable (Methods, 2017).")


def _ligated_construct() -> Construct:
    return Construct([*_copy_segments(), Segment("3-prime dA", "A"), _a_core(), _a_tail()],
                     name="one-end-ligated R-ChIP DNA")


def _ligated_scene() -> Scene:
    copy = _copy_segments()
    sc = Scene()
    sc.strand("copy", [*copy, Segment("3-prime dA", "A"), _a_core(), _a_tail()],
              label="ligated copy strand", mod5="invddT")
    sc.anneal("template", complement_segments(copy[1:]), to="copy",
              pair=("R-loop DNA-derived insert'", "R-loop DNA-derived insert"),
              label="captured template DNA")
    sc.anneal("oligo B", [Segment("B 5-prime tail", OLIGO_B[:8]),
                            Segment("B paired core", OLIGO_B[8:20]),
                            Segment("B 3-prime T", OLIGO_B[20:])],
              to="copy", pair=("B paired core", "oligo A paired core"),
              label="adapter oligo B", unpaired=("B 5-prime tail",))
    sc.junction("copy", "3-prime dA", "oligo A paired core", "adapter ligation")
    sc.note("template", "B 3′-T / template 5′ nick: sealing needs a 5′-phosphorylated template end (not stated)")
    return sc


def _rchip_library() -> Construct:
    return Construct([
        *_copy_segments(),
        Segment("dA-tail-derived A", "A"),
        Segment("adapter oligo A", OLIGO_A, "r1"),
        Segment("5-nt flexible barcode (barcode primer)", "N" * 5, "cbc",
                placeholder=True, feature=SAMPLE_BARCODE),
        Segment("P5 (barcode primer)", revcomp(il.P5), "p5"),
    ], name="strand-specific R-ChIP library")


def _oligo(name: str, seq: str) -> Construct:
    return Construct([Segment(name, seq)], name=name)


DEFINITION = Definition("R-ChIP — in-vivo R-loop capture with inactive RNASEH1",
        "10.1016/j.molcel.2017.10.008", "Chen et al., Molecular Cell (2017)",
        "V5-tagged RNASEH1 D210N binds R-loops in fixed cells; a tailed N9 primer and one-ended adapter ligation produce a strand-specific DNA-component library.",
        ("tailed N9 end", "R-loop DNA-derived insert"),
        "Single-end 40-nt Read 1 (Illumina Read 1 primer, HiSeq 2500) starts at the ligated adapter and reads the captured DNA strand's sequence.",
        "01_r_chip.html",
        caveat="The final library is reconstructed from the published N9 primer, adapter strands and PCR primers (2019 protocol, Table 1). How the 5-nt barcode between P5 and Read 1 was read is not stated in the inspected sources.")

FINAL_LIBRARY = _rchip_library()
READ1 = replace(sp.TRUSEQ["R1"], note=(
    "2019 protocol, ‘Sequencing and basic data analysis’: “Sequence the R-ChIP library with Illumina read 1 sequencing "
    "primer on the HiSeq 2500”; the primer sequence itself is Illumina's."))
SEQ_PRIMERS = (READ1,)
READ_LENGTHS = {"Read 1": 40}
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Detailed chemistry was checked against the authors’ 2019 Nature Protocols protocol (PMC6604627).")
SPEC.CAVEAT = DEFINITION.caveat
SPEC.READ_LENGTHS = READ_LENGTHS
SPEC.oligos = lambda: (
    f"<b>Tailed N9 primer</b> — 5′-/invddt/{N9_PRIMER}-3′",
    f"<b>Adapter oligo A</b> — 5′-/Phos/{OLIGO_A}-3′",
    f"<b>Adapter oligo B</b> — 5′-{OLIGO_B}-3′",
    f"<b>PCR primer</b> — 5′-{PCR_PRIMER}-3′",
    f"<b>Barcoded PCR primer</b> — 5′-{BARCODE_PRIMER}-3′ (NNNNN: flexible barcode)",
)


def adapter_geometry():
    return rchip_adapter_geometry()


def _hybrid_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed graph: 2017 Methods (“Strand-specific R-ChIP Library Construction”) with
    reagents and conditions from the authors' 2019 protocol (PMC6604627, Steps 61–74 onward and its sequencing section)."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    cap = RLoopCapture('V5-tagged catalytically inactive RNASEH1 D210N', 'DNA', in_vivo=True)
    chromatin = graph.input('rloop_chromatin', 'chromatin containing R-loops', MaterialKind.BIOLOGICAL, state=cap.rloop(), primary=True)
    formaldehyde = graph.input('formaldehyde', '1% formaldehyde, 10–15 min RT; glycine quench', MaterialKind.CONDITION)
    fixed = graph.product('fixed_chromatin', 'cross-linked R-loop chromatin', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('crosslink', 'Cross-link cells', StepKind.CROSS_LINKING, ((chromatin, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    sonication = graph.input('sonication', 'nuclear lysis and sonication to 250–600 bp', MaterialKind.CONDITION)
    sheared = graph.product('sheared_chromatin', 'sheared fixed chromatin', MaterialKind.MOLECULAR, state=cap.rloop())
    graph.react('shear', 'Shear cross-linked chromatin', StepKind.FRAGMENTATION, ((fixed, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (sheared,))
    v5_antibody = graph.input('v5_antibody', 'anti-V5 antibody magnetic beads (overnight, 4 °C)', MaterialKind.SELECTION)
    inactive_rnaseh = graph.input('inactive_rnaseh1', 'V5–RNASEH1 D210N R-loop binder', MaterialKind.SELECTION)
    captured = graph.product('captured_rloop_dna', 'captured R-loop chromatin', MaterialKind.MOLECULAR, state=cap.captured_component())
    rejected = graph.product('uncaptured_chromatin', 'uncaptured chromatin (washed off with TSE I–III, TE)', MaterialKind.WASTE)
    graph.react('immunoprecipitate', 'Immunoprecipitate V5–RNASEH1-bound R-loops', StepKind.AFFINITY_CAPTURE, ((sheared, InputRole.TEMPLATE), (inactive_rnaseh, InputRole.SELECTION), (v5_antibody, InputRole.SELECTION)), (captured, rejected))
    graph.discard(rejected, 'not retained by anti-V5 immunoprecipitation; washed away')
    reverse = graph.input('decrosslink_condition', 'SDS elution; overnight 65 °C de-cross-linking', MaterialKind.CONDITION)
    rnase_pk = graph.input('rnase_a_proteinase_k', 'RNase A then proteinase K', MaterialKind.CATALYST)
    purified = graph.product('purified_hybrid_dna', 'purified precipitated hybrid fragments', MaterialKind.MOLECULAR, state=cap.captured_component())
    digested = graph.product('protein_rna_digest', 'digested protein and RNase A-sensitive RNA, extracted', MaterialKind.WASTE)
    graph.react('decrosslink_purify', 'Elute, de-cross-link, digest and phenol-extract', StepKind.DEPROTEINIZATION, ((captured, InputRole.TEMPLATE), (reverse, InputRole.CONDITION), (rnase_pk, InputRole.CATALYST)), (purified, digested))
    graph.discard(digested, 'removed by phenol and phenol:chloroform:isoamyl alcohol extraction and ethanol precipitation')
    reviewed = ['rloop_chromatin', 'formaldehyde', 'crosslink', 'sonication', 'shear', 'v5_antibody', 'inactive_rnaseh1', 'immunoprecipitate', 'decrosslink_condition', 'rnase_a_proteinase_k', 'decrosslink_purify']
    denaturing = graph.input('denaturing_condition', '95 °C 5 min with N9 primer, then 25 °C 5 min', MaterialKind.CONDITION)
    single = graph.product('single_stranded_captured_dna', 'single-stranded captured DNA with annealed N9 primer', MaterialKind.MOLECULAR, state=cap.captured_component())
    n9 = graph.input('tailed_n9_primer', 'tailed N9 random primer (5′-/invddt/)', MaterialKind.OLIGO, state=_oligo('tailed N9 primer', N9_PRIMER))
    graph.react('denature', 'Denature captured DNA and anneal the tailed N9 primer', StepKind.HYBRIDIZATION, ((purified, InputRole.TEMPLATE), (n9, InputRole.STRUCTURAL), (denaturing, InputRole.CONDITION)), (single,))
    polymerase = graph.input('copy_polymerase', 'phi29 DNA polymerase (30 °C 20 min; 65 °C 10 min)', MaterialKind.CATALYST)
    primed_scene = Scene()
    primed_scene.strand('copy', _copy_segments(), label='tailed-N9 copy strand', mod5='invddT')
    primed_scene.anneal('template', complement_segments(_copy_segments()[1:]), to='copy', pair=("R-loop DNA-derived insert'", 'R-loop DNA-derived insert'), label='captured template DNA')
    primed = graph.product('n9_copied_dna', 'tailed-N9 copied dsDNA', MaterialKind.MOLECULAR, state=primed_scene)
    graph.react('random_primer_copy', 'Extend the tailed N9 primer on the captured DNA', StepKind.EXTENSION, ((single, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (primed,))
    tailing = graph.input('a_tailing_enzyme', 'Klenow fragment (3′→5′ exo−), dATP, 37 °C 30 min', MaterialKind.CATALYST)
    tailed = graph.product('one_end_da_tailed', 'one-end dA-tailed copied DNA', MaterialKind.MOLECULAR, state=_one_end_tailed_scene())
    graph.react('da_tail', 'Add one 3-prime dA to the unobstructed end', StepKind.TAILING, ((primed, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input('illumina_adapter', 'annealed oligo A/B adapter (Table 1)', MaterialKind.ADAPTER, state=_adapter_scene())
    ligase = graph.input('dna_ligase', 'T4 DNA ligase, 2 h RT', MaterialKind.CATALYST)
    tagged = graph.product('one_end_ligated_dna', 'one-end-adapter-ligated DNA', MaterialKind.MOLECULAR, state=_ligated_scene())
    graph.react('adapter_ligation', 'Ligate the adapter to the dA-tailed end only', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,), geometry=rchip_adapter_geometry())
    reviewed.extend(('denaturing_condition', 'denature', 'tailed_n9_primer', 'copy_polymerase', 'random_primer_copy', 'a_tailing_enzyme', 'da_tail', 'illumina_adapter', 'dna_ligase', 'adapter_ligation'))
    forward = graph.input('library_forward_primer', 'PCR primer (Table 1)', MaterialKind.OLIGO, state=_oligo('PCR primer', PCR_PRIMER))
    reverse_p = graph.input('library_reverse_primer', 'barcode primer (Table 1)', MaterialKind.OLIGO, state=_oligo('barcode primer', BARCODE_PRIMER))
    pcr_polymerase = graph.input('pcr_polymerase', 'Phusion High-Fidelity DNA polymerase', MaterialKind.CATALYST)
    amplified = graph.product('amplified_library', 'amplified R-ChIP library', MaterialKind.MOLECULAR, state=final)
    graph.react('library_pcr', 'PCR-amplify (test 14/16 cycles; 16–18 in the 2017 Methods)', StepKind.PCR, ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse_p, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.exact(_ligated_construct(), PCRPrimer('PCR primer', PCR_PRIMER), PCRPrimer('barcode primer', BARCODE_PRIMER),
                    statement='The PCR primer repeats the N9-primer tail; the barcode primer’s 3′ TruSeq Read 1 end anneals to ligated oligo A plus the dA (2019 protocol, Table 1).'))
    gel = graph.input('gel_isolation', '2% agarose gel isolation (see size-range note)', MaterialKind.SELECTION)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final)
    off_size = graph.product('off_size_products', 'products outside the excised size range', MaterialKind.WASTE)
    graph.react('size_select', 'Gel-isolate the library size range', StepKind.SELECTION, ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, off_size))
    graph.discard(off_size, 'left in the agarose gel')
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS, acquired_roles=('Read 1',), read_lengths=READ_LENGTHS, read_uses={'Read 1': 'include'}))
    graph.annotate(SourceDiscrepancy(
        'size_select', 'Excised library size range',
        (SourceClaim('Chen et al. 2017 Mol Cell, STAR Methods “Strand-specific R-ChIP Library Construction”', '16–18 PCR cycles; products of 130–350 bp gel-isolated'),
         SourceClaim('Chen et al. 2019 Nat Protoc (PMC6604627), PCR and library recovery steps', 'test PCR at 14 and 16 cycles (18 if weak); gel slice of 140–350 bp excised')),
        'The schematic does not choose a range or cycle number; both authors’ values are shown.',
        'The two documents are separate author publications describing the same library step.',
        'Which range applied to a given deposited dataset is not resolved by these passages.'))
    reviewed.extend(('library_forward_primer', 'library_reverse_primer', 'pcr_polymerase', 'library_pcr', 'amplified_library', 'gel_isolation', 'size_select', 'sequencing_library'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _hybrid_graph(d, final)
SPEC.protocol_graph = protocol_graph
