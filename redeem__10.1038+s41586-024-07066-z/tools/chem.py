"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from batch_ngs import nextera_library
from batch_ngs import truseq_library
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_specs import role
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("ReDeeM — deep mitochondrial mutation profiling", "10.1038/s41586-024-07066-z", "Weng et al., Nature (2024)", "A whole-cell 10x Multiome modification retains mtDNA, then tiled hybrid capture produces a deeply sequenced mtDNA branch matched to RNA and ATAC.", "truseq", ("cell barcode", "mtDNA fragment"), "The mtDNA, ATAC and RNA libraries carry matched cell barcodes; mtDNA fragment endpoints support consensus error correction.", "multi")


from protocol_specs import workflow, two_branch


def _source_workflow():
    return workflow(two_branch("fixed permeabilized whole cell", "nuclear chromatin", "retained mtDNA"), (
        ("Run 10x Multiome while retaining mitochondrial DNA",
         two_branch("matched cell barcode", "RNA + ATAC libraries", "mtDNA fragments"),
         "All branches retain a common droplet cell identity."),
        ("Hybrid-capture mtDNA in four staggered probe pools",
         molecule("captured mtDNA branch", ("cell barcode", "mtDNA fragment", "capture probe")),
         "Published Supplementary Data 1 lists the probe sequences; fragment endpoints act as endogenous molecular identifiers."),
        ("Sequence mtDNA, ATAC and RNA as three matched libraries",
         two_branch("same cell barcode", "mtDNA library", "ATAC + RNA libraries"),
         "The mtDNA branch receives deep error-corrected sequencing."),
    ))

SPEC = get('redeem', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "whole cells retaining nuclear chromatin, RNA and mtDNA",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("whole cell", ("nuclear chromatin", "RNA", "mtDNA")),
                        primary=True)
    droplets = graph.input("multiome_droplets", "10x Multiome droplet/barcoding condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "matched cell-barcoded RNA, ATAC and mtDNA products",
                             MaterialKind.MOLECULAR,
                             state=molecule("matched cell products", (
                                 "cell barcode", "RNA products", "ATAC fragments", "mtDNA fragments")))
    failed = graph.product("failed_partitions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("multiome", "Run 10x Multiome while retaining mitochondrial DNA",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable shared cell barcode")
    mt = graph.product("mt_branch", "cell-barcoded mitochondrial-DNA branch",
                       MaterialKind.MOLECULAR, state=barcoded.state)
    atac = graph.product("atac_branch", "matched ATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    rna = graph.product("rna_branch", "matched RNA branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    graph.react("split_branches", "Separate mtDNA, ATAC and RNA branches",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (mt, atac, rna),
                changes_molecule=False)
    probes = graph.input("capture_probes", "four staggered mtDNA capture-probe pools (sequences in Supplementary Data 1)",
                         MaterialKind.SELECTION)
    captured = graph.product("captured_mt", "hybrid-captured cell-barcoded mtDNA",
                             MaterialKind.MOLECULAR, state=mt.state)
    rejected = graph.product("rejected_mt", "non-mitochondrial capture material",
                             MaterialKind.WASTE)
    graph.react("capture_mt", "Hybrid-capture mtDNA in four staggered probe pools",
                StepKind.AFFINITY_CAPTURE,
                ((mt, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not belong to the selected mitochondrial branch")
    reviewed = ["cells", "multiome_droplets", "multiome", "split_branches",
                "capture_probes", "capture_mt"]
    finish_pcr(graph, "mt", captured, final, primers, "mtDNA branch: " + d.sequencing,
                reviewed, label="ReDeeM mtDNA")
    atac_final, atac_primers = nextera_library(
        [role("cell barcode"), role("accessible-DNA insert")], "ReDeeM ATAC library")
    rna_final, rna_primers = truseq_library(
        [role("cell barcode"), role("UMI"), role("transcript cDNA")], "ReDeeM RNA library")
    finish_pcr(graph, "atac", atac, atac_final, atac_primers,
                "matched ATAC branch: " + d.sequencing, reviewed, label="ReDeeM ATAC")
    finish_pcr(graph, "rna", rna, rna_final, rna_primers,
                "matched RNA branch: " + d.sequencing, reviewed, label="ReDeeM RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
