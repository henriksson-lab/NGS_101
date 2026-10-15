"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("MALBAC-DT — high-detectability single-cell transcriptomics", "10.1101/2019.12.31.892190", "Chen et al., bioRxiv (2020; later PNAS 2022)", "Barcoded oligo(dT) reverse transcription is followed by MALBAC quasi-linear random amplification, pooling, tagmentation and indexed sequencing.", "nextera", ("cell barcode", "UMI", "transcript cDNA"), "Read 1 reports the cell barcode and UMI; Read 2 reports transcript sequence after TruePrep tagmentation.", "RNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("single-cell poly(A) RNA", ("transcripts",)), (
        ("Reverse-transcribe with a cell barcode, UMI and poly(dT)",
         molecule("barcoded first-strand cDNA", ("cell barcode", "UMI", "cDNA"), joined=True),
         "The well barcode and UMI are installed before pooling."),
        ("Perform MALBAC quasi-linear random-primer amplification",
         molecule("MALBAC-DT amplicons", ("common handle", "barcoded cDNA", "common handle")),
         "Looping suppresses repeated copying of early amplicons."),
        ("Pool, tagment and add sequencing indices",
         molecule("indexed library", ("left Nextera end", "cell barcode", "UMI", "cDNA", "right Nextera end"), joined=True),
         "Read 1 reports barcode and UMI; Read 2 reports transcript sequence."),
    ))

SPEC = get('malbac-dt', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "single-cell poly(A) RNA", MaterialKind.BIOLOGICAL,
                      state=molecule("single-cell poly(A) RNA", ("transcripts",)), primary=True)
    rt_primer = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "cell-barcoded first-strand cDNA", MaterialKind.MOLECULAR,
                         state=molecule("barcoded first-strand cDNA", (
                             "cell barcode", "UMI", "transcript cDNA")))
    graph.react("reverse_transcribe", "Reverse-transcribe with barcode, UMI and oligo(dT)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    random_primer = graph.input("malbac_primer", "MALBAC common-handle random primer",
                                MaterialKind.OLIGO)
    strand_polymerase = graph.input("malbac_polymerase", "strand-displacing MALBAC polymerase",
                                    MaterialKind.CATALYST)
    amplified = graph.product("amplified", "quasi-linearly amplified MALBAC-DT cDNA",
                              MaterialKind.MOLECULAR,
                              state=molecule("MALBAC-DT amplicons", (
                                  "common handle", "barcoded cDNA", "common handle")))
    graph.react("malbac", "Perform MALBAC quasi-linear random-primer amplification",
                StepKind.EXTENSION,
                ((cdna, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (strand_polymerase, InputRole.CATALYST)), (amplified,))
    transposome = graph.input("transposome", "TruePrep loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "TruePrep transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented MALBAC-DT products", MaterialKind.MOLECULAR,
                           state=amplified.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Pool and tagment the amplified cDNA", StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable barcoded transcript insert")
    reviewed = ["rna", "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "malbac_primer", "malbac_polymerase", "malbac", "transposome", "tn5", "tagment"]
    finish_pcr(graph, "library", tagged, final, primers, d.sequencing, reviewed,
                label="MALBAC-DT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
