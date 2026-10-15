from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, custom_library, single, token
from rna_cell_assays import feature

from rna_cell_assays import chain

def lineage_workflow() -> Workflow:
    return chain("integrated scGESTALT transgene", single("lineage cassette", (
        token("DsRed transcript"), token("CRISPR target array"), token("poly(A)"))), (
        ("Induce Cas9 and guide pools at successive developmental times",
         single("edited lineage cassette", (token("DsRed transcript"),
                                              token("cumulative target-array indels"), token("poly(A)"))),
         "Repair-derived insertions and deletions form a permanent, cumulative lineage record."),
        ("Capture cassette RNA and transcriptome by inDrops",
         single("cell-indexed cDNA", (token("transcript or lineage cDNA"),
                                      token("cell-barcode part 1", tag="cbc", length=8),
                                      token("cell-barcode part 2", tag="cbc", length=8),
                                      token("UMI", tag="umi", length=6))),
         "The inDrops cell identifier is the join key between expression and lineage libraries."),
        ("IVT-amplify, then split fragmented and unfragmented RNA",
         branches("same indexed amplified RNA", (
             ("transcriptome", (token("fragmented amplified RNA"),)),
             ("lineage", (token("unfragmented lineage RNA"),)))),
         "Fragmented RNA makes the expression library; unfragmented RNA preserves the full target array."),
        ("Nested-PCR the target array and add matched Illumina indexes",
         single("scGESTALT lineage library", (token("lineage barcode array"),
                                               token("cell-barcode part 1", tag="cbc", length=8),
                                               token("cell-barcode part 2", tag="cbc", length=8),
                                               token("UMI", tag="umi", length=6), token("library index", tag="cbc", length=8))),
         "The lineage and transcriptome libraries use the same multiplexing index so cells can be matched."),
    ))

DEFINITION = Definition("scGESTALT — single-cell lineage recording", "10.1038/s41596-018-0058-x",
        "Raj &amp; Gagnon, <i>Nature Protocols</i> (2018)", "01_scgestalt.html",
        "Cumulative CRISPR indels in an expressed target array record ancestry; inDrops captures that lineage transcript alongside the cell transcriptome and supplies their shared cell identifier.",
        "The lineage run uses Read 1 = 250 cycles for the edited array, Index 1 = 8 cycles for cell-barcode part 1, Index 2 = 8 cycles for library index, and Read 2 = 14 cycles for cell-barcode part 2 plus 6-nt UMI.")

d = DEFINITION
final = custom_library('scGESTALT lineage library', (token('lineage barcode array', length=44), token('cell-barcode part 1', tag='cbc', length=8, molecular_feature=feature('cell_barcode_1', 'cell_barcode', 'combinatorial', group='cell_barcode', part='part 1')), token('cell-barcode part 2', tag='cbc', length=8, molecular_feature=feature('cell_barcode_2', 'cell_barcode', 'combinatorial', group='cell_barcode', part='part 2')), token('UMI', tag='umi', length=6, molecular_feature=feature('umi', 'umi', 'random')), token('library index', tag='cbc', length=8, molecular_feature=feature('sample_index', 'sample_index', 'fixed'))))
primers = ()
reads = None
OLIGOS = ()
workflow_fn = lineage_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def scgestalt_graph(definition, final: Construct) -> ProtocolGraph:
    """scGESTALT graph with matched transcriptome and lineage-library outcomes."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input(
        "transgenic_cell", "cell with integrated expressed lineage cassette",
        MaterialKind.BIOLOGICAL,
        state=single("lineage cassette and transcriptome", (
            token("DsRed transcript"), token("CRISPR target array"),
            token("poly(A)"), token("cellular transcriptome"))), primary=True)
    guides = graph.input("guide_pools", "successive CRISPR guide pools",
                         MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "inducible Cas9", MaterialKind.CATALYST)
    edited_cell = graph.product(
        "edited_cell", "cell with cumulatively edited lineage cassette",
        MaterialKind.MOLECULAR,
        state=single("edited lineage cassette and transcriptome", (
            token("DsRed transcript"), token("cumulative target-array indels"),
            token("poly(A)"), token("cellular transcriptome"))))
    graph.react("edit_lineage", "Induce cumulative lineage-cassette editing",
                StepKind.GENOME_EDITING,
                ((cell, InputRole.TEMPLATE), (guides, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (edited_cell,))
    host_rnap = graph.input("host_rna_polymerase", "cellular RNA polymerase",
                            MaterialKind.CATALYST)
    cell_rna = graph.product(
        "cellular_rna", "polyadenylated edited-cassette RNA plus cellular mRNA",
        MaterialKind.MOLECULAR,
        state=single("polyadenylated cellular RNA", (
            token("edited lineage transcript"), token("cellular transcripts"))))
    graph.react("transcribe_cassette", "Transcribe and polyadenylate the edited cassette",
                StepKind.TRANSCRIPTION,
                ((edited_cell, InputRole.TEMPLATE),
                 (host_rnap, InputRole.CATALYST)), (cell_rna,))

    indrops_primer = graph.input("indrops_capture_primer",
                                 "inDrops cell-barcode/UMI oligo(dT) primer",
                                 MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    first_strand = graph.product(
        "barcoded_first_strand", "cell-indexed lineage and transcriptome cDNA",
        MaterialKind.MOLECULAR,
        state=single("cell-indexed cDNA", (
            token("transcript or lineage cDNA"),
            token("cell-barcode part 1", tag="cbc", length=8),
            token("cell-barcode part 2", tag="cbc", length=8),
            token("UMI", tag="umi", length=6))))
    graph.react("reverse_transcription", "Capture lineage RNA and transcriptome by inDrops",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell_rna, InputRole.TEMPLATE),
                 (indrops_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first_strand,))
    second_primer = graph.input("second_strand_primer", "inDrops second-strand primer",
                                MaterialKind.OLIGO, inferred=True)
    second_polymerase = graph.input("second_strand_polymerase", "second-strand polymerase",
                                    MaterialKind.CATALYST)
    duplex = graph.product("barcoded_cdna_duplex", "barcoded second-strand cDNA",
                           MaterialKind.MOLECULAR, state=first_strand.state)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((first_strand, InputRole.TEMPLATE),
                 (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (duplex,))
    t7 = graph.input("t7_rna_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    amplified_rna = graph.product("amplified_rna", "IVT-amplified barcoded RNA",
                                  MaterialKind.MOLECULAR, state=duplex.state)
    graph.react("ivt", "Linearly amplify barcoded cDNA by IVT", StepKind.TRANSCRIPTION,
                ((duplex, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (amplified_rna,))
    transcriptome_rna = graph.product("transcriptome_rna", "RNA for transcriptome library",
                                      MaterialKind.MOLECULAR,
                                      state=single("transcriptome branch", (
                                          token("amplified transcriptome RNA"),)))
    lineage_rna = graph.product("lineage_rna", "unfragmented lineage-cassette RNA",
                                MaterialKind.MOLECULAR,
                                state=single("lineage branch", (
                                    token("unfragmented lineage RNA"),)))
    graph.react("split_ivt", "Split amplified RNA into transcriptome and lineage branches",
                StepKind.SPLIT, ((amplified_rna, InputRole.TEMPLATE),),
                (transcriptome_rna, lineage_rna))

    fragmentation = graph.input("rna_fragmentation_reagent", "RNA fragmentation reagent",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_transcriptome_rna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR, state=transcriptome_rna.state)
    graph.react("fragment_transcriptome", "Fragment the transcriptome branch",
                StepKind.FRAGMENTATION,
                ((transcriptome_rna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))
    rna_adapter = graph.input("transcriptome_rna_adapter", "inDrops RNA-library adapter",
                              MaterialKind.ADAPTER, inferred=True)
    rna_ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated_rna = graph.product("adapter_ligated_transcriptome_rna",
                                "adapter-ligated transcriptome RNA",
                                MaterialKind.MOLECULAR, state=fragmented.state,
                                inferred=True)
    graph.react("transcriptome_adapter_ligation", "Ligate the transcriptome RNA adapter",
                StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE),
                 (rna_adapter, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (ligated_rna,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The accessible scGESTALT protocol establishes the inDrops branch but not one exact adapter-end sequence."))
    transcriptome_rt_primer = graph.input("transcriptome_rt_primer",
                                          "transcriptome RT primer",
                                          MaterialKind.OLIGO, inferred=True)
    expression_cdna = graph.product("expression_cdna", "transcriptome-library cDNA",
                                    MaterialKind.MOLECULAR, state=ligated_rna.state,
                                    inferred=True)
    graph.react("transcriptome_rt", "Reverse-transcribe the fragmented transcriptome RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated_rna, InputRole.TEMPLATE),
                 (transcriptome_rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (expression_cdna,))
    pcr_polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    expression_forward = graph.input("expression_forward_primer",
                                     "inDrops expression PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "inDrops expression PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_final = custom_library("scGESTALT matched transcriptome library", (
        token("cell barcode", tag="cbc", length=16), token("UMI", tag="umi", length=6),
        token("transcript cDNA", length=44)))
    expression_library = graph.product(
        "expression_library", "matched scGESTALT transcriptome library",
        MaterialKind.LIBRARY, state=expression_final, inferred=True)
    graph.react("expression_pcr", "PCR the transcriptome branch",
                StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The accessible protocol states the matched inDrops run but not one assembled primer pair."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The transcriptome branch follows the matched inDrops library implementation."))

    lineage_rt_primer = graph.input("lineage_random_hexamer", "lineage random-hexamer RT primer",
                                    MaterialKind.OLIGO)
    lineage_cdna = graph.product("lineage_cdna", "full target-array lineage cDNA",
                                 MaterialKind.MOLECULAR, state=lineage_rna.state)
    graph.react("lineage_rt", "Reverse-transcribe unfragmented lineage RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lineage_rna, InputRole.TEMPLATE),
                 (lineage_rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (lineage_cdna,))
    nested1_forward = graph.input("lineage_nested1_forward", "lineage nested-PCR-1 forward primer",
                                  MaterialKind.OLIGO)
    nested1_reverse = graph.input("lineage_nested1_reverse", "lineage nested-PCR-1 reverse primer",
                                  MaterialKind.OLIGO)
    nested = graph.product("lineage_nested_amplicon", "first nested lineage amplicon",
                           MaterialKind.MOLECULAR, state=lineage_cdna.state)
    graph.react("lineage_nested_pcr1", "First nested PCR of the lineage array",
                StepKind.PCR,
                ((lineage_cdna, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The page preserves published role geometry; exact assembled primer sites remain unavailable."))
    nested2_forward = graph.input("lineage_nested2_forward", "indexed lineage PCR forward primer",
                                  MaterialKind.OLIGO)
    nested2_reverse = graph.input("lineage_nested2_reverse", "indexed lineage PCR reverse primer",
                                  MaterialKind.OLIGO)
    lineage_library = graph.product("lineage_library", "scGESTALT lineage library",
                                    MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("lineage_nested_pcr2", "Second nested PCR adds matched Illumina indexes",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE),
                 (nested2_forward, InputRole.STRUCTURAL),
                 (nested2_reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (lineage_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Exact inDrops primer bases are role-labelled in the accessible protocol."))
    graph.sequence(lineage_library, SequencingPlan(
        limitation=definition.endpoint,
        read_lengths={"Read 1": 250, "Index 1": 8, "Index 2": 8, "Read 2": 14}))
    graph.declare_reviewed_scope((
        "transgenic_cell", "guide_pools", "cas9", "edit_lineage",
        "host_rna_polymerase", "transcribe_cassette", "indrops_capture_primer",
        "reverse_transcriptase", "reverse_transcription", "second_strand_primer",
        "second_strand_polymerase", "second_strand", "t7_rna_polymerase", "ivt",
        "split_ivt", "rna_fragmentation_reagent", "fragment_transcriptome",
        "transcriptome_rna_adapter", "rna_ligase", "transcriptome_adapter_ligation",
        "transcriptome_rt_primer", "transcriptome_rt", "pcr_polymerase",
        "expression_forward_primer", "expression_reverse_primer", "expression_pcr",
        "expression_library", "lineage_random_hexamer", "lineage_rt",
        "lineage_nested1_forward", "lineage_nested1_reverse", "lineage_nested_pcr1",
        "lineage_nested2_forward", "lineage_nested2_reverse", "lineage_nested_pcr2",
        "lineage_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'scgestalt'
    return scgestalt_graph(d, final)
SPEC.protocol_graph = protocol_graph
