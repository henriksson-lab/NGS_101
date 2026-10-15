"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from circular import CircularizationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from rna_library_operations import molecule
from stlfr_chem import barcode_transfer_geometry as stlfr_barcode_transfer_geometry

DEFINITION = Definition("stLFR — single-tube long-fragment reads", "10.1101/gr.245126.118", "Wang et al., Genome Research (2019)", "Transposon-tagged subfragments from one long DNA molecule are captured on one bead and receive a shared linked-read barcode before DNBSEQ.", "dnb", ("linked-read barcode", "genomic subfragment"), "The cited detailed protocol publishes the BGI Read 1, barcode and Read 2 sequencing primers and states paired-end 100-base genomic reads; barcode demultiplexing consumes a 42-base barcode/common-adapter string. The exact DNB circle-closing adapter or splint geometry is not published.", "DNA", "The paper-cited Protocol Exchange workbook publishes the common bead, bridge, transposon, PCR, branch-ligation, BGI and Illumina oligos, plus all 1,536 members of each of the three 10-nt barcode sets. The page shows the barcode-family architecture; the external workbook remains the authoritative full whitelist.")


from protocol_specs import workflow
from stlfr_chem import oligos as stlfr_oligos


def _source_workflow():
    return workflow(molecule("20–300-kb DNA", ("long genomic molecule",)), (
        ("Insert transposons carrying capture sequences along each long molecule",
         molecule("tagged long DNA", ("subfragment", "capture transposon", "subfragment")),
         "Transposition creates many addressable sites without losing long-range provenance."),
        ("Capture one long molecule on one barcoded bead",
         molecule("bead-bound molecule", ("bead barcode", "captured subfragments")),
         "Millions of bead barcodes provide practically nonredundant co-barcoding."),
        ("Ligate the bead barcode to captured subfragments and make DNBs",
         molecule("co-barcoded library", ("linked-read barcode", "short genomic subfragment"), joined=True,
                  junctions=(("linked-read barcode", "short genomic subfragment", "barcode ligation"),)),
         "Reads sharing a barcode are assigned back to the same long molecule."),
    ))

SPEC = get('stlfr', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=stlfr_oligos,
           sequencing_unavailable='The BGI sequencing-primer sequences are published; only the DNB circle-closing adapter or splint geometry remains unavailable.')


def _local_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "20–300-kb genomic DNA molecules", MaterialKind.BIOLOGICAL,
                      state=molecule("long genomic molecule", ("long genomic molecule",)),
                      primary=True)
    transposome = graph.input("transposome", "capture-sequence-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "long molecules carrying addressable transposons",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagged long DNA", (
                               "subfragment", "capture transposon", "subfragment")))
    waste = graph.product("tn5_waste", "non-library transposition products", MaterialKind.WASTE)
    graph.react("tagment", "Insert capture transposons along each long molecule",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an addressable genomic subfragment")
    bead = graph.input("barcoded_bead", "single barcoded capture bead", MaterialKind.SELECTION)
    captured = graph.product("captured", "one long molecule captured per barcoded bead",
                             MaterialKind.MOLECULAR,
                             state=molecule("bead-bound molecule", (
                                 "bead barcode", "captured subfragments")))
    rejected = graph.product("uncaptured", "uncaptured or multiply assigned fragments",
                             MaterialKind.WASTE)
    graph.react("capture", "Capture one long molecule on one barcoded bead",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE), (bead, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not have an unambiguous long-molecule bead assignment")
    barcode_adapter = graph.input("barcode_adapter", "bead-barcode ligating oligo",
                                  MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "barcode DNA ligase", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "co-barcoded genomic subfragments",
                             MaterialKind.MOLECULAR,
                             state=molecule("co-barcoded library precursor", (
                                 "linked-read barcode", "short genomic subfragment"),
                                 junctions=(("linked-read barcode", "short genomic subfragment",
                                             "barcode ligation"),)))
    graph.react("ligate_barcode", "Ligate the bead barcode to captured subfragments",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (barcode_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (barcoded,),
                geometry=stlfr_barcode_transfer_geometry())
    circle_ligase = graph.input("circle_ligase", "DNB library circularization ligase",
                                MaterialKind.CATALYST)
    library = graph.product("library", "stLFR DNBSEQ library", MaterialKind.LIBRARY,
                            state=final)
    graph.react("circularize", "Circularize the co-barcoded products for DNB formation",
                StepKind.CIRCULARIZATION,
                ((barcoded, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)),
                (library,), geometry=CircularizationGeometry.unavailable(
                    "The stLFR source identifies DNB library circularization but does not disclose the circle-closing adapter/splint geometry."))
    graph.sequence(library, SequencingPlan(limitation=d.sequencing))
    reviewed = ["dna", "transposome", "tn5", "tagment", "barcoded_bead", "capture",
                "barcode_adapter", "ligase", "ligate_barcode", "circle_ligase",
                "circularize", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY)

SPEC.protocol_graph = protocol_graph
