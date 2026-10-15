from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, single, token
from rna_cell_assays import custom_library

from rna_cell_assays import chain

def run_on_workflow() -> Workflow:
    """Nascent-RNA run-on with the label at the extended 3' end."""
    substrate = "isolated nuclei"
    labelled = "BrU run-on end"
    selection = "anti-BrdU immunopurification is repeated three times"
    return chain(substrate, single(substrate, (token("engaged RNA polymerase"),
                                               token("nascent RNA"))), (
        ("Allow engaged polymerases to run on; prevent new initiation",
         single("extended nascent RNA", (token("pre-existing nascent RNA"),
                                          token(labelled, tag="umi", length=8))),
         "The nucleotide analogue is incorporated only at the RNA 3′ end extended in vitro."),
        ("Fragment RNA and select analogue-containing molecules",
         single("selected nascent RNA", (token("nascent RNA fragment"),
                                          token(labelled, tag="umi", length=8))), selection),
        ("Install the RNA-library adapter and reverse-transcribe",
         single("nascent-RNA cDNA", (token("library handle"),
                                      token("nascent-RNA-derived cDNA"))),
         "The genomic end of the read reports the position and orientation of engaged polymerase."),
        ("PCR to a sequencing-compatible library",
         single("run-on library", (token("left sequencing arm"),
                                    token("nascent-RNA-derived cDNA"),
                                    token("right sequencing arm"))),
         "The defining source does not expose one universal modern primer set for this library."),
    ))

DEFINITION = Definition("GRO-seq — global run-on sequencing", "10.1126/science.1162228",
        "Core, Waterfall &amp; Lis, <i>Science</i> (2008)", "01_gro-seq.html",
        "A nuclear run-on incorporates BrUTP into short extensions made by already-engaged polymerases; repeated anti-BrdU purification selects those nascent RNAs for strand-specific sequencing.",
        "The original paper used 33-base single-end Illumina reads from the 5′ end of fragmented, BrU-labelled RNA; its historical library does not expose a current canonical primer set.")

d = DEFINITION
final = custom_library(d.title + ' library', (token('nascent-RNA-derived cDNA', length=44),))
primers = ()
reads = None
OLIGOS = ()
workflow_fn = run_on_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def run_on_graph(*, chromatin: bool, definition, final: Construct) -> ProtocolGraph:
    """Typed run-on graph shared by GRO-seq and ChRO-seq."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    substrate_label = "isolated chromatin" if chromatin else "isolated nuclei"
    label = "biotin-NTP" if chromatin else "BrUTP"
    substrate = graph.input(
        "run_on_substrate", substrate_label, MaterialKind.BIOLOGICAL,
        state=single(substrate_label, (token("engaged RNA polymerase"),
                                       token("nascent RNA"))), primary=True)
    polymerase = graph.input("engaged_polymerase", "already-engaged RNA polymerase",
                             MaterialKind.CATALYST)
    analogue = graph.input("labelled_ntp", label + " run-on nucleotide",
                           MaterialKind.MOLECULAR)
    labelled = graph.product(
        "labelled_nascent_rna", "run-on-labelled nascent RNA",
        MaterialKind.MOLECULAR,
        state=single("extended nascent RNA", (
            token("pre-existing nascent RNA"),
            token(label + " run-on end", tag="umi", length=8))))
    graph.react("run_on", "Extend already-engaged nascent transcripts",
                StepKind.TRANSCRIPTION,
                ((substrate, InputRole.TEMPLATE),
                 (polymerase, InputRole.CATALYST),
                 (analogue, InputRole.STRUCTURAL)), (labelled,))

    fragmentation = graph.input("rna_fragmentation_reagent",
                                "RNA fragmentation reagent", MaterialKind.CONDITION)
    fragments = graph.product("labelled_rna_fragments", "fragmented labelled nascent RNA",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment_rna", "Fragment run-on-labelled RNA", StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    selector = graph.input(
        "label_selector",
        "streptavidin" if chromatin else "anti-BrdU immunopurification reagent",
        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_nascent_rna", "selected analogue-containing nascent RNA",
        MaterialKind.MOLECULAR,
        state=single("selected nascent RNA", (
            token("nascent RNA fragment"),
            token(label + " run-on end", tag="umi", length=8))))
    rejected = graph.product("unlabelled_rna", "unselected RNA", MaterialKind.WASTE)
    graph.react("select_labelled_rna", "Select analogue-containing RNA",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE),
                 (selector, InputRole.SELECTION)), (selected, rejected))
    graph.discard(rejected, "lacks the run-on label required for affinity selection")

    adapter = graph.input("rna_adapter", "run-on RNA-library adapter",
                          MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated_rna", "adapter-ligated nascent RNA",
                            MaterialKind.MOLECULAR, state=selected.state,
                            inferred=True)
    graph.react("adapter_ligation", "Install the RNA-library adapter",
                StepKind.LIGATION,
                ((selected, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining historical implementation does not expose one universal modern RNA-adapter end."))
    rt_primer = graph.input("rt_primer", "adapter-complementary RT primer",
                            MaterialKind.OLIGO, inferred=True)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product(
        "nascent_cdna", "nascent-RNA-derived cDNA", MaterialKind.MOLECULAR,
        state=single("nascent-RNA cDNA", (
            token("library handle"), token("nascent-RNA-derived cDNA"))),
        inferred=True)
    graph.react("reverse_transcription", "Reverse-transcribe selected nascent RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cdna,))
    forward = graph.input("pcr_forward_primer", "library PCR forward primer",
                          MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("pcr_reverse_primer", "library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True)
    pcr_enzyme = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", definition.title + " library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "PCR to a sequencing-compatible library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (pcr_enzyme, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The defining source does not expose one universal modern primer pair."))
    graph.sequence(library, SequencingPlan(limitation=definition.endpoint))
    graph.declare_reviewed_scope((
        "run_on_substrate", "engaged_polymerase", "labelled_ntp", "run_on",
        "rna_fragmentation_reagent", "fragment_rna", "label_selector",
        "select_labelled_rna", "rna_adapter", "rna_ligase", "adapter_ligation",
        "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "pcr_forward_primer", "pcr_reverse_primer", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'gro-seq'
    return run_on_graph(chromatin=False, definition=d, final=final)
SPEC.protocol_graph = protocol_graph
