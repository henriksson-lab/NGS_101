"""Local source and reaction model for truseq-stranded-mrna."""
from chemdraw import Construct, Scene, Segment, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from endprep import cohesive_ligation_geometry, dA_tailed_scene, ligation_product_scene, repair_and_dA_tail
from illumina_vendor_workflows import Definition, TRUSEQ_Y_END, _finish_pcr, _single, _truseq_final, _workflow, truseq_final_rows
from types import SimpleNamespace


def _stranded_mrna_workflow(final: Construct) -> Workflow:
    selected = _single("RNA", [Segment("fragmented poly(A) RNA", "R" * 36,
                                       placeholder=True)], label="fragmented poly(A) RNA")
    first = _single("first strand", [Segment("antisense first-strand cDNA", "N" * 36,
                                              placeholder=True)], label="first-strand cDNA")
    marked = Scene.duplex([
        Segment("antisense first-strand cDNA", "N" * 36, placeholder=True)
    ], label="dUTP-marked cDNA duplex")
    marked.mark("bottom", "antisense first-strand cDNA'", "second strand contains dUTP")
    end_prep = repair_and_dA_tail(Construct([
        Segment("stranded cDNA insert", "N" * 36, placeholder=True)]))
    tailed = dA_tailed_scene(end_prep, "dA-tailed stranded cDNA")
    w = _workflow("total RNA", _single("RNA", [Segment("total RNA", "R" * 40,
                                                        placeholder=True)], label="total RNA"), (
        ("Select poly(A) RNA twice and fragment at elevated temperature",
         "fragmented poly(A) RNA", selected,
         "Poly-T magnetic beads enrich polyadenylated RNA; divalent cations fragment it."),
        ("Random-prime first-strand reverse transcription",
         "antisense first-strand cDNA", first,
         "Actinomycin D suppresses DNA-dependent synthesis during first-strand synthesis."),
        ("RNase H / DNA polymerase I second-strand synthesis with dUTP",
         "dUTP-marked duplex", marked,
         "Replacing dTTP with dUTP marks the second strand for exclusion during enrichment."),
        ("End-repair and add one 3-prime dA", "dA-tailed cDNA", tailed,
         "A-tailed cDNA accepts the complementary T-overhang adapters."),
    ))
    w.ligate_adapter(
        "Ligate TruSeq adapters and selectively PCR-enrich",
        cohesive_ligation_geometry(end_prep, TRUSEQ_Y_END,
            statement="The public TruSeq adapter has a 3′-T end complementary to dA-tailed cDNA."),
        ligation_product_scene(final, "antisense cDNA insert", label="stranded library",
                               adapter_group="TruSeq adapter").rows(),
        name="stranded sequencing library",
        note="The enrichment polymerase does not copy past dUTP, so Read 1 maps antisense and Read 2 sense.")
    return w

def _stranded_mrna_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "total RNA", MaterialKind.BIOLOGICAL,
                      state=_single("RNA", [Segment("total RNA", "R" * 40,
                                                   placeholder=True)], label="total RNA"),
                      primary=True)
    beads = graph.input("polyt_beads", "poly-T RNA-selection beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "poly(A)-selected RNA", MaterialKind.MOLECULAR,
                             state=_single("RNA", [Segment("poly(A) RNA", "R" * 36,
                                                         placeholder=True)], label="poly(A) RNA"))
    rejected = graph.product("non_polya", "non-polyadenylated RNA", MaterialKind.WASTE)
    graph.react("select_polya", "Select poly(A) RNA twice", StepKind.AFFINITY_CAPTURE,
                ((rna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "is depleted by the poly-T selection")
    fragmentation = graph.input("fragmentation_condition", "elevated-temperature divalent-cation fragmentation condition",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented", "fragmented poly(A) RNA",
                               MaterialKind.MOLECULAR, state=selected.state)
    graph.react("fragment", "Fragment selected poly(A) RNA", StepKind.FRAGMENTATION,
                ((selected, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragmented,))
    random_primer = graph.input("rt_primer", "random first-strand RT primer",
                                MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase with actinomycin D",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "antisense first-strand cDNA",
                          MaterialKind.MOLECULAR,
                          state=_single("first strand", [Segment("antisense first-strand cDNA",
                              "N" * 36, placeholder=True)], label="first-strand cDNA"))
    graph.react("reverse_transcribe", "Random-prime first-strand reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragmented, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    second_polymerase = graph.input("second_strand_polymerase", "RNase H / DNA polymerase I mix",
                                    MaterialKind.CATALYST)
    marked = graph.product("marked_duplex", "dUTP-marked cDNA duplex",
                           MaterialKind.MOLECULAR, state=first.state)
    graph.react("second_strand", "Synthesize the dUTP-marked second strand",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (marked,))
    repair = graph.input("repair_enzyme", "cDNA end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt stranded-cDNA duplexes",
                             MaterialKind.MOLECULAR, state=marked.state)
    graph.react("end_repair", "End-repair the cDNA duplex", StepKind.END_REPAIR,
                ((marked, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme",
                          MaterialKind.CATALYST)
    end_prep = repair_and_dA_tail(Construct([
        Segment("stranded cDNA insert", "N" * 36, placeholder=True)]))
    tailed = graph.product("tailed", "dA-tailed stranded cDNA", MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(end_prep, "dA-tailed stranded cDNA"))
    graph.react("tail", "Add one 3-prime dA to the cDNA ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("adapter", "TruSeq 3-prime-T Y adapters", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "TruSeq adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated directional cDNA",
                            MaterialKind.MOLECULAR, state=tailed.state)
    geometry = cohesive_ligation_geometry(
        end_prep, TRUSEQ_Y_END,
        statement="The public TruSeq adapter has a 3′-T end complementary to dA-tailed cDNA.")
    graph.react("ligate", "Ligate TruSeq adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=geometry)
    reviewed = ["rna", "polyt_beads", "select_polya", "fragmentation_condition", "fragment",
                "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "second_strand_primer", "second_strand_polymerase", "second_strand",
                "repair_enzyme", "end_repair", "tailing_enzyme", "tail", "adapter",
                "ligase", "ligate"]
    _finish_pcr(graph, "library", adapted, library, primers, d.sequencing_intro,
                reviewed, label="TruSeq stranded mRNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "TruSeq Stranded mRNA", "01_truseq_stranded_mrna.html",
        'Illumina <a href="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-stranded-mrna-workflow/truseq-stranded-mrna-workflow-reference-1000000040498-00.pdf">TruSeq Stranded mRNA Reference Guide 1000000040498 v00</a>.',
        "Poly(A) RNA is bead-selected and fragmented, random-prime copied, marked with dUTP in the second strand, adapter-ligated and selectively PCR-enriched.",
        "The guide supports single- and dual-index variants. The schematic uses the dual-index HS architecture and does not assign a particular index identity.",
        "Directional TruSeq cDNA library: Read 1 derives from antisense cDNA and Read 2 from the sense transcript strand.",
        "Standard TruSeq primers sequence the antisense-derived insert in Read 1, the sense-derived mate in Read 2 and both sample indexes.",
        "truseq", "antisense cDNA insert", _stranded_mrna_workflow)


def protocol_graph():
    return _stranded_mrna_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)

FINAL_LIBRARY, SEQ_PRIMERS = _truseq_final(DEFINITION.insert_name, DEFINITION.title + " library")
d = DEFINITION
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
    CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
    workflow=lambda: _stranded_mrna_workflow(FINAL_LIBRARY),
    final_rows=lambda: truseq_final_rows(FINAL_LIBRARY),
    protocol_graph=protocol_graph,
)
globals().update(vars(SPEC))
