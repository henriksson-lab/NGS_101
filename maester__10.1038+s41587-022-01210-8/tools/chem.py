"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from pcr import PCRGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("MAESTER — mitochondrial variant enrichment from scRNA-seq", "10.1038/s41587-022-01210-8", "Miller et al., Nature Biotechnology (2022)", "Multiplex PCR enriches mitochondrial transcripts from barcoded 3′ scRNA-seq cDNA while preserving cell barcodes and UMIs.", "truseq", ("cell barcode", "UMI", "mitochondrial target"), "Read 1 retains the parent scRNA-seq cell barcode and UMI; the opposite read reports the enriched mitochondrial target.", "RNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("barcoded 3′ scRNA-seq cDNA", ("cell barcode", "UMI", "cDNA")), (
        ("Multiplex-PCR mitochondrial transcript targets",
         molecule("targeted amplicon", ("cell barcode", "UMI", "mitochondrial target"), joined=True),
         "Target enrichment preserves the original cell barcode and UMI."),
        ("Pool target reactions and add Illumina indices",
         molecule("MAESTER library", ("sample index", "cell barcode", "UMI", "mtRNA target"), joined=True),
         "The enriched library is sequenced separately from the parent expression library."),
    ))

SPEC = get('maester', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cdna = graph.input("cdna", "barcoded 3-prime single-cell RNA-seq cDNA",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("barcoded scRNA-seq cDNA", ("cell barcode", "UMI", "cDNA")),
                       primary=True)
    target_forward = graph.input("target_forward", "multiplex mitochondrial target primers",
                                 MaterialKind.OLIGO)
    target_reverse = graph.input("target_reverse", "parent-library reverse primer",
                                 MaterialKind.OLIGO)
    polymerase = graph.input("target_polymerase", "target PCR polymerase",
                             MaterialKind.CATALYST)
    enriched = graph.product("enriched", "cell-barcode/UMI-preserving mitochondrial amplicons",
                             MaterialKind.MOLECULAR,
                             state=molecule("targeted amplicon", (
                                 "cell barcode", "UMI", "mitochondrial target")))
    graph.react("target_pcr", "Multiplex-PCR mitochondrial transcript targets",
                StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (target_forward, InputRole.STRUCTURAL),
                 (target_reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (enriched,), geometry=PCRGeometry.unavailable(
                    enriched.state, reason="The multiplex target panel, rather than one universal amplicon, defines the PCR geometry."))
    reviewed = ["cdna", "target_forward", "target_reverse", "target_polymerase", "target_pcr"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="MAESTER mitochondrial-target")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
