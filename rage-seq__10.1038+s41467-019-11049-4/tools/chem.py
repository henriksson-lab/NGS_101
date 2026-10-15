"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import truseq_library
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from protocol_specs import role
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("RAGE-Seq — repertoire and gene expression by sequencing", "10.1038/s41467-019-11049-4", "Singh et al., Nature Communications (2019)", "Hybrid capture enriches full-length immune-receptor cDNAs from barcoded droplet libraries for long-read sequencing alongside short-read expression.", "nanopore", ("motor adapter", "cell barcode", "full-length receptor cDNA"), "The long read retains the cell barcode and full receptor isoform; the matched short-read library supplies gene expression. The paper prints the TSO and R1 primers used before and after both capture rounds; Oxford Nanopore kit adapters remain proprietary.", "RNA", "The original Methods and supplement were rechecked. Exact TSO and R1 amplification primers are published. The 644-exon Roche–NimbleGen capture panel is specified by target genes rather than probe sequences, and the SQK-LSK108/SQK-LSK308 adapter sequences are not disclosed.")


from protocol_specs import workflow
from rage_seq_chem import oligos as rage_seq_oligos


def _source_workflow():
    return workflow(molecule("droplet scRNA-seq cDNA", ("cell barcode", "UMI", "full-length cDNA")), (
        ("Hybrid-capture immunoglobulin and T-cell-receptor cDNA",
         molecule("captured receptor cDNA", ("cell barcode", "UMI", "full-length V(D)J transcript")),
         "Biotinylated probes enrich receptor molecules while retaining cellular identifiers."),
        ("Attach long-read adapters and sequence the captured branch",
         molecule("RAGE-Seq long-read molecule", ("Nanopore adapter", "cell barcode", "full-length receptor cDNA"), joined=True,
                  junctions=(("Nanopore adapter", "cell barcode", "adapter ligation"),)),
         "Short-read gene expression and long-read receptor data join through the barcode."),
    ))

SPEC = get('rage-seq', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=rage_seq_oligos,
           sequencing_mechanism="nanopore_motor",
           motor_adapter_name="SQK-LSK108 / SQK-LSK308 motor adapter")


def _local_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cdna = graph.input("cdna", "droplet scRNA-seq cDNA carrying cell barcode and UMI",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("droplet scRNA-seq cDNA", (
                           "cell barcode", "UMI", "full-length cDNA")), primary=True)
    expression = graph.product("expression_branch", "short-read expression branch",
                               MaterialKind.MOLECULAR, state=cdna.state)
    receptor = graph.product("receptor_branch", "immune-receptor enrichment branch",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("split", "Split expression and immune-receptor branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),), (expression, receptor),
                changes_molecule=False)
    probes = graph.input("capture_probes", "biotinylated immunoglobulin/TCR probes",
                         MaterialKind.SELECTION)
    captured = graph.product("captured", "captured full-length receptor cDNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured receptor cDNA", (
                                 "cell barcode", "UMI", "full-length receptor transcript")))
    rejected = graph.product("rejected", "non-receptor cDNA in the capture branch",
                             MaterialKind.WASTE)
    graph.react("capture", "Hybrid-capture immunoglobulin and TCR cDNA",
                StepKind.AFFINITY_CAPTURE,
                ((receptor, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not contain a selected immune-receptor transcript")
    motor = graph.input("motor_adapter", "Nanopore motor sequencing adapter",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "Nanopore adapter ligase", MaterialKind.CATALYST)
    long_library = graph.product("long_library", "RAGE-Seq long-read receptor library",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("ligate_motor", "Attach long-read adapters to captured receptor cDNA",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (long_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Nanopore adapter end chemistry is kit-revision-specific; the source establishes its molecular role."))
    graph.sequence(long_library, SequencingPlan(limitation=d.sequencing))
    expression_final, expression_primers = truseq_library(
        [role("cell barcode"), role("UMI"), role("transcript cDNA")],
        "RAGE-Seq matched expression library")
    reviewed = ["cdna", "split", "capture_probes", "capture", "motor_adapter",
                "adapter_ligase", "ligate_motor", "long_library"]
    finish_pcr(graph, "expression", expression, expression_final, expression_primers,
                "matched short-read expression branch", reviewed,
                label="RAGE-Seq expression")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY)

SPEC.protocol_graph = protocol_graph
