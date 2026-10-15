from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, single, sp, token
from rna_cell_assays import nextera_library

from rna_cell_assays import chain

def find_seq_workflow() -> Workflow:
    return chain("heterogeneous cell suspension", single("one cell", (token("genome"),
                                                                        token("poly(A) RNA"))), (
        ("Co-encapsulate one cell, molten allyl agarose, and acrydite oligo(dT)",
         single("hydrogel bead precursor", (token("cell genome"), token("cell RNA"),
                                             token("covalent oligo(dT) matrix"))),
         "APS/TEMED links acrydite oligo(dT) to allyl agarose."),
        ("Lyse strongly and reverse-transcribe mRNA on the hydrogel",
         single("genome-capture bead", (token("retained genome"), token("matrix-bound cDNA"))),
         "Proteinase K/LiDS exposes DNA while the covalent capture matrix retains the transcriptome."),
        ("Re-encapsulate beads, run marker-specific TaqMan PCR, and sort droplets",
         single("fluorescent selected bead", (token("retained genome"), token("matrix-bound cDNA"),
                                               token("TaqMan marker signal"))),
         "Marker PCR is a selection measurement; it is not the transcriptome library."),
        ("Recover selected cDNA and build Smart-seq2 or Smart-seq3 libraries",
         single("FIND-seq library", (token("Nextera left arm"), token("selected-cell cDNA"),
                                      token("Nextera right arm"))),
         "The 2024 protocol supports bulk bins or individually sorted beads and UMI-aware Smart-seq3."),
    ))

DEFINITION = Definition("FIND-seq — nucleic-acid cytometry for rare cells", "10.1038/s41596-024-01021-y",
        "Shin et al., <i>Nature Protocols</i> (2024)", "01_find-seq.html",
        "A hydrogel bead retains one cell's genome and oligo(dT)-captured transcriptome while marker-specific TaqMan PCR drives fluorescent droplet sorting before whole-transcriptome library preparation.",
        "Recovered selected-cell cDNA is prepared by Smart-seq2 or UMI-aware Smart-seq3 and sequenced through its validated Nextera library ends.")

d = DEFINITION
final, primers = nextera_library([token('selected-cell transcript cDNA', length=44)], 'FIND-seq transcriptome library')
reads = None
OLIGOS = (
            "<b>Smart-seq2 bead capture</b> — 5′-acrydite-TTTTTTAAGCAGTGGTATCAACGCAGAGTAC(T)₃₀VN-3′",
            "<b>Smart-seq2 template switch</b> — 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′",
            "<b>Smart-seq3 bead capture</b> — 5′-acrydite-TTTTTTACGAGCATCAGCAGCATACGA(T)₃₀V-3′",
            "<b>Smart-seq3 template switch</b> — 5′-biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′; N₈ is the UMI",
        )
workflow_fn = find_seq_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def find_seq_graph(definition, final: Construct,
                   primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """FIND-seq graph from hydrogel partition through selected-cell library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cells = graph.input("cell_suspension", "heterogeneous cell suspension",
                        MaterialKind.BIOLOGICAL, state=find_seq_workflow().initial, primary=True)
    matrix = graph.input("allyl_agarose", "molten allyl agarose matrix",
                         MaterialKind.CONDITION)
    capture_oligo = graph.input("acrydite_oligodt", "acrydite oligo(dT)",
                                MaterialKind.OLIGO)
    bead = graph.product("cell_hydrogel", "one-cell hydrogel bead precursor",
                         MaterialKind.MOLECULAR,
                         state=single("hydrogel bead precursor", (
                             token("cell genome"), token("cell RNA"),
                             token("covalent oligo(dT) matrix"))))
    empty = graph.product("empty_partition_material", "empty or multiplet partitions",
                          MaterialKind.WASTE)
    graph.react("encapsulate", "Encapsulate cells in oligo(dT) hydrogels",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (matrix, InputRole.CONDITION),
                 (capture_oligo, InputRole.STRUCTURAL)), (bead, empty))
    graph.discard(empty, "partitions without one usable cell are not selected for sequencing")
    proteinase = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    lids = graph.input("lids", "LiDS lysis condition", MaterialKind.CONDITION)
    lysed = graph.product("lysed_bead", "lysed genome-retaining hydrogel bead",
                          MaterialKind.MOLECULAR, state=bead.state)
    graph.react("lyse", "Lyse cells strongly inside the hydrogel",
                StepKind.DEPROTEINIZATION,
                ((bead, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST),
                 (lids, InputRole.CONDITION)), (lysed,))
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    copied = graph.product(
        "genome_capture_bead", "genome-retaining bead with matrix-bound cDNA",
        MaterialKind.MOLECULAR,
        state=single("genome-capture bead", (
            token("retained genome"), token("matrix-bound cDNA"))))
    graph.react("reverse_transcription", "Reverse-transcribe mRNA on the hydrogel",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lysed, InputRole.TEMPLATE),
                 (capture_oligo, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    taqman = graph.input("taqman_assay", "marker-specific TaqMan assay",
                         MaterialKind.SELECTION)
    selected = graph.product(
        "selected_bead", "fluorescence-selected genome-capture bead",
        MaterialKind.MOLECULAR,
        state=single("fluorescent selected bead", (
            token("retained genome"), token("matrix-bound cDNA"),
            token("TaqMan marker signal"))))
    rejected = graph.product("rejected_beads", "marker-negative beads",
                             MaterialKind.WASTE)
    graph.react("marker_sort", "Run marker-specific TaqMan PCR and sort droplets",
                StepKind.AFFINITY_CAPTURE,
                ((copied, InputRole.TEMPLATE), (taqman, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "marker-negative droplets are excluded by fluorescence sorting")
    recovery = graph.input("cdna_recovery_condition", "selected-bead cDNA recovery",
                           MaterialKind.CONDITION)
    recovered = graph.product("recovered_cdna", "selected-cell cDNA",
                              MaterialKind.MOLECULAR,
                              state=single("selected-cell cDNA", (token("full-length cDNA"),)))
    graph.react("recover_cdna", "Recover selected-cell cDNA", StepKind.EXTRACTION,
                ((selected, InputRole.TEMPLATE),
                 (recovery, InputRole.CONDITION)), (recovered,))
    preamp_forward = graph.input("smart_forward_primer", "Smart-seq amplification forward primer",
                                 MaterialKind.OLIGO, inferred=True)
    preamp_reverse = graph.input("smart_reverse_primer", "Smart-seq amplification reverse primer",
                                 MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("library_polymerase", "library amplification polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "Smart-seq-amplified selected-cell cDNA",
                              MaterialKind.MOLECULAR, state=recovered.state,
                              inferred=True)
    graph.react("smartseq_pcr", "Amplify selected cDNA by Smart-seq2 or Smart-seq3",
                StepKind.PCR,
                ((recovered, InputRole.TEMPLATE),
                 (preamp_forward, InputRole.STRUCTURAL),
                 (preamp_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The protocol permits Smart-seq2 or UMI-aware Smart-seq3."))
    loaded = graph.input("loaded_nextera_adapters", "loaded Nextera adapters",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=amplified.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment selected-cell cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (loaded, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    gap_fill = graph.input("gap_fill_polymerase", "Nextera gap-fill polymerase",
                           MaterialKind.CATALYST)
    filled = graph.product("gap_filled_cdna", "gap-filled tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (gap_fill, InputRole.CATALYST)), (filled,))
    index_forward = graph.input("index_forward_primer", "Nextera indexed forward primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse_primer", "Nextera indexed reverse primer",
                                MaterialKind.OLIGO)
    library = graph.product("sequencing_library", "FIND-seq transcriptome library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Index and amplify the FIND-seq library", StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Index choice varies; the final Nextera endpoint is validated."))
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.declare_reviewed_scope((
        "cell_suspension", "allyl_agarose", "acrydite_oligodt", "encapsulate",
        "proteinase_k", "lids", "lyse", "reverse_transcriptase",
        "reverse_transcription", "taqman_assay", "marker_sort",
        "cdna_recovery_condition", "recover_cdna", "smart_forward_primer",
        "smart_reverse_primer", "library_polymerase", "smartseq_pcr",
        "loaded_nextera_adapters", "tn5", "tagment", "gap_fill_polymerase",
        "fill_tagmentation_gaps", "index_forward_primer",
        "index_reverse_primer", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'find-seq'
    return find_seq_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
