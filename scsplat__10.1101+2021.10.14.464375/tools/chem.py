"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_graph import SourceClaim, SourceDiscrepancy
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("scSPLAT — scalable single-cell WGBS", "10.1101/2021.10.14.464375", "Raine et al., bioRxiv (2021; Scientific Reports 2022)", "Bisulfite-converted single-cell DNA receives a P5/cell-barcode random-primer copy before pooling and bulk splinted ligation of the opposite adapter.", "truseq", ("cell barcode", "bisulfite-converted insert"), "The first six bases of Read 1 are the inline cell barcode; paired genomic reads report the converted insert.", "DNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("one cell", ("genomic DNA",)), (
        ("Bisulfite-convert DNA separately in each well",
         molecule("single-stranded converted DNA", ("bisulfite-converted fragments",)),
         "Conversion precedes adapter installation."),
        ("Random-prime second-strand synthesis with P5 and an inline cell barcode",
         molecule("barcoded strand", ("P5", "cell barcode", "random-primer copy"), joined=True),
         "Up to 32 individually barcoded wells can now be pooled."),
        ("Pool cells and install the other adapter by splinted ligation",
         molecule("scSPLAT library", ("P5", "cell barcode", "converted insert", "splint-ligated P7"), joined=True,
                  junctions=(("converted insert", "splint-ligated P7", "splinted ligation"),)),
         "The splint makes adapter tagging a bulk reaction after cell indexing."),
    ))

SPEC = get('scsplat', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "single-cell genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("one cell", ("genomic DNA",)), primary=True)
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product("converted", "single-stranded bisulfite-converted DNA",
                              MaterialKind.MOLECULAR,
                              state=molecule("converted DNA", ("bisulfite-converted fragments",)))
    graph.react("convert", "Bisulfite-convert DNA in each well",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
    random_primer = graph.input("random_primer", "P5/cell-barcode random primer",
                                MaterialKind.OLIGO)
    polymerase = graph.input("strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "P5/cell-barcode second-strand copy",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded strand", (
                                 "P5", "cell barcode", "random-primer copy")))
    graph.react("second_strand", "Random-prime second-strand synthesis with P5 and cell barcode",
                StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (barcoded,))
    adapter = graph.input("p7_adapter", "splinted P7 adapter", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "splinted-adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "pooled splint-ligated single-cell templates",
                            MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("ligate", "Pool wells and install P7 by splinted ligation",
                StepKind.LIGATION,
                ((barcoded, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source identifies splinted P7 installation; this role-level page does not assert the panel-specific splint bases."))
    reviewed = ["dna", "bisulfite", "convert", "random_primer", "strand_polymerase",
                "second_strand", "p7_adapter", "ligase", "ligate"]
    finish_pcr(graph, "library", adapted, final, primers, d.sequencing, reviewed,
                label="scSPLAT")
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "p7_adapter", "Published SPLAT adapter-name conflict",
        (SourceClaim("Raine et al., Scientific Reports (2022), Methods, adapter annealing",
                     "The adapter is described as annealed from ss1a and ss2b."),
         SourceClaim("Same Methods, printed SPLAT-adapter oligo list",
                     "The second oligo is named ss1b: 5′Phos/AGATCGGAAGAGCACACGTC/3′AmMo.")),
        "The graph records a splinted P7 adapter with unavailable exact end geometry; it does not choose whether ss1b and ss2b are the same reagent.",
        "The Methods print the second oligo and its modifications but use two names for it.",
        "Supplementary Table S1 prints ss1b and no ss2b entry; the intended annealing-partner name remains uncertain."))
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
