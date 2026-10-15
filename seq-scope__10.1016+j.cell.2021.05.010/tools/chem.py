"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("Seq-Scope — flow-cell spatial transcriptomics", "10.1016/j.cell.2021.05.010", "Cho et al., Cell (2021)", "An Illumina flow cell is first decoded as a dense HDMI barcode array and then reused to capture spatially resolved tissue transcripts.", "truseq", ("spatial HDMI", "UMI", "transcript cDNA"), "The first run uses a published custom Read1-DraI primer to decode HDMI coordinates. The second paired-end run uses the completed TruSeq shell: one read reports HDMI/UMI and the other reports captured cDNA.", "RNA", "The open Methods print all primary and alternative seed oligos, both custom first-run primers, the random second-strand primer and both PCR rounds with phosphorothioates. These were previously missed.")


from protocol_specs import workflow
from seqscope_chem import oligos as seqscope_oligos


def _source_workflow():
    return workflow(molecule("Illumina flow cell", ("random HDMI oligos",)), (
        ("Cluster and sequence the HDMI array in situ",
         molecule("decoded capture surface", ("mapped HDMI cluster", "poly(dT) capture end")),
         "The first sequencing run maps each spatial barcode to an imaging coordinate."),
        ("Place tissue, capture mRNA and reverse-transcribe on the flow cell",
         molecule("surface cDNA", ("spatial HDMI", "UMI", "transcript cDNA"), joined=True),
         "Each cDNA inherits the barcode of the cluster beneath it."),
        ("Release, amplify and index the spatial transcript library",
         molecule("Seq-Scope library", ("sample index", "HDMI", "UMI", "transcript"), joined=True),
         "Sequencing reconnects transcript identity with the decoded surface coordinate."),
    ))

SPEC = get('seq-scope', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=seqscope_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    surface = graph.input("surface", "Illumina flow cell bearing random HDMI capture oligos",
                          MaterialKind.BIOLOGICAL,
                          state=molecule("flow-cell array", ("random HDMI oligos",)), primary=True)
    decoder = graph.input("decoder", "in-situ HDMI sequencing primers",
                          MaterialKind.OLIGO)
    decoded = graph.product("decoded", "coordinate-decoded HDMI capture surface",
                            MaterialKind.MOLECULAR,
                            state=molecule("decoded capture surface", (
                                "mapped HDMI cluster", "poly(dT) capture end")))
    graph.react("decode", "Cluster and sequence the HDMI array in situ",
                StepKind.HYBRIDIZATION,
                ((surface, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (decoded,))
    tissue_rna = graph.input("tissue_rna", "poly(A) RNA released from overlaid tissue",
                             MaterialKind.MOLECULAR)
    captured = graph.product("captured", "RNA hybridized to decoded surface clusters",
                             MaterialKind.MOLECULAR, state=decoded.state)
    graph.react("capture", "Capture tissue mRNA on the decoded flow cell",
                StepKind.HYBRIDIZATION,
                ((decoded, InputRole.TEMPLATE), (tissue_rna, InputRole.STRUCTURAL)),
                (captured,))
    rt_primer = graph.input("rt_primer", "surface-bound poly(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "surface cDNA carrying HDMI and UMI",
                         MaterialKind.MOLECULAR,
                         state=molecule("surface cDNA", ("spatial HDMI", "UMI", "transcript cDNA")))
    graph.react("reverse_transcribe", "Reverse-transcribe captured RNA on the surface",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed = ["surface", "decoder", "decode", "tissue_rna", "capture", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe"]
    finish_pcr(graph, "library", cdna, final, primers, d.sequencing, reviewed,
                label="Seq-Scope spatial transcript")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
