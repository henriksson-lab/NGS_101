"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("SNuBar-ATAC — oligo sample multiplexing for scATAC", "10.1016/j.molcel.2021.09.026", "Wang et al., Molecular Cell (2021)", "A sample-barcode oligo anneals to the universal transposome tail during tagmentation before samples are pooled into 10x scATAC droplets.", "nextera", ("cell barcode", "sample barcode", "UMI", "accessible DNA"), "The standard 10x ATAC reads recover genomic insert and cell barcode; the SBO-derived read reports sample barcode and UMI.", "DNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("sample-specific nuclei", ("chromatin",)), (
        ("Tagment each sample with the standard 10x ATAC transposome",
         molecule("tagmented chromatin", ("genomic fragment", "transposome universal tail")),
         "The sample barcode is not carried by Tn5 itself."),
        ("Anneal a SNuBar oligo to the transposome tail",
         molecule("sample-tagged fragment", ("genomic fragment", "sample barcode", "UMI", "PCR handle"), joined=True),
         "The SBO transfers sample identity before samples are pooled."),
        ("Pool samples and run the 10x single-cell ATAC workflow",
         molecule("SNuBar-ATAC library", ("10x cell barcode", "sample barcode", "UMI", "ATAC fragment"), joined=True),
         "Cell and sample identities occupy distinct molecular fields."),
    ))

SPEC = get('snubar-atac', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "sample-separated nuclei containing accessible chromatin",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("sample-specific nuclei", ("accessible chromatin",)),
                         primary=True)
    transposome = graph.input("transposome", "standard 10x ATAC loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented chromatin carrying the universal tail",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagmented chromatin", (
                               "genomic fragment", "transposome universal tail")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment each sample with the 10x ATAC transposome",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an accessible-DNA insert with the universal tail")
    snubar = graph.input("snubar", "sample-barcode/UMI SNuBar oligo",
                         MaterialKind.OLIGO)
    sample_tagged = graph.product("sample_tagged", "SNuBar-tagged ATAC fragments",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("sample-tagged fragment", (
                                      "genomic fragment", "sample barcode", "UMI", "PCR handle")))
    graph.react("anneal_snubar", "Anneal the SNuBar oligo to the transposome tail",
                StepKind.HYBRIDIZATION,
                ((tagged, InputRole.TEMPLATE), (snubar, InputRole.STRUCTURAL)),
                (sample_tagged,))
    droplets = graph.input("droplets", "10x single-cell ATAC droplet condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell- and sample-barcoded ATAC molecules",
                             MaterialKind.MOLECULAR,
                             state=molecule("SNuBar-ATAC products", (
                                 "10x cell barcode", "sample barcode", "UMI", "ATAC fragment")))
    collisions = graph.product("collisions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Pool samples and partition nuclei into 10x droplets",
                StepKind.PARTITION,
                ((sample_tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, collisions), changes_molecule=False)
    graph.discard(collisions, "does not retain an interpretable cell barcode")
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "snubar", "anneal_snubar",
                "droplets", "partition"]
    finish_pcr(graph, "library", barcoded, final, primers, d.sequencing, reviewed,
                label="SNuBar-ATAC")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
