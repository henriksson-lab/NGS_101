"""MGI/DNBSEQ ss-circle and DNA-nanoball preparation model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, circle_rows, feature
from circular import CircularizationGeometry, circularize_ssdna, rolling_circle
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "DNBSEQ library circularization and DNA-nanoball preparation"
NOTES = "01_dnbseq.html"
SOURCE = 'Commercial protocol: MGI <a href="https://en.mgi-tech.com/Download/download_file_new/id/441">MGIEasy Fast FS Library Prep Set User Manual 2.0</a> and DNBSEQ DNB-preparation documentation.'
SUMMARY = "An indexed double-stranded library is denatured, one strand is splint-ligated into a single-stranded circle, remaining linear DNA is digested, and rolling-circle replication compacts tandem copies into a DNA nanoball for patterned-array sequencing."
CAVEAT = "MGI adapter, splint and sequencing-primer bases are proprietary in the cited current kit documents. Their roles and topology are shown without fabricated sequences."
SEQ_PRIMERS = ()
FINAL_LIBRARY = None
SEQUENCING_ENDING = "A DNBSEQ sequencing primer anneals to the repeated adapter site within each DNA nanoball. Because rolling-circle replication creates many tandem copies of the same ssDNA circle, every repeat presents the same primer site and insert boundary to cPAS sequencing. The current primer bases are proprietary, so binding geometry—not an invented sequence—is shown."

LINEAR = Construct([seg("adapter A / primer site", "X" * 16, "r1", placeholder=True,
                        inferred=True),
                    seg("insert", "N" * 30, placeholder=True),
                    seg("adapter B / barcode", "Y" * 16, "cbc", placeholder=True,
                        inferred=True, feature=feature("sample_barcode_region", "sample_index", "unknown",
                                                       note="barcode sub-boundary is proprietary"))], name="denatured DNBSEQ library strand")
CIRCLE = circularize_ssdna(LINEAR, five_prime_phosphate=True)
RCA = rolling_circle(CIRCLE, "X" * 16, copies=3)

def sections():
    upstream=AdapterLigationGeometry.unavailable("The current MGI kit manuals describe upstream adapter ligation but keep the adapter strands and ligation-facing end geometry proprietary.")
    circular=AdapterLigationGeometry.unavailable("The DNB-preparation manual describes splint-assisted circle closure but does not disclose the splint bases or the aligned adapter-end sequences.")
    return [
        ("Denature the indexed double-stranded library",
         upstream.rows()+[Row(chunks=[("indexed dsDNA → heat/alkali → one adapter-bearing ssDNA strand", None, False)])],
         "Upstream fragmentation, repair, adapter ligation and PCR produce the platform library; DNB preparation begins by making it single-stranded."),
        ("Splint-ligate the strand into a closed circle",
         circular.rows()+circle_rows(CIRCLE.linear, closure_label="splint ligation"),
         "A complementary splint aligns the two adapter ends for ligation. The closure joins the final adapter segment back to the first; exonuclease removes unclosed linear molecules."),
        ("Copy the circle into a compact DNA nanoball",
         [Row(chunks=[("circle → φ29 rolling-circle replication", None, False)]),
          Row(chunks=[("[adapter—insert—barcode]×[adapter—insert—barcode]×[adapter—insert—barcode]", "r1", False)]),
          Row(chunks=[("                  ↓ self-compaction into one DNB", None, False)])],
         "Rolling-circle replication produces a single concatemer containing tandem complementary copies of the original circle; it is compacted and loaded onto a patterned array."),
        ("Prime repeated adapter sites for cPAS sequencing",
         [Row(chunks=[("…[primer site]—insert—barcode—[primer site]—insert—barcode…", "r1", False)]),
          Row(chunks=[("   sequencing primer →              sequencing primer →", "r2", False)])],
         "The same undisclosed sequencing-primer site occurs once per RCA repeat. Probe-anchor synthesis interrogates the immobilized nanoball."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    ds_library = graph.input("indexed_library", "adapter-bearing indexed dsDNA library",
                             MaterialKind.MOLECULAR, primary=True)
    denaturant = graph.input("denaturant", "heat or alkaline denaturation condition",
                             MaterialKind.CONDITION)
    strand = graph.product("library_strand", LINEAR.name, MaterialKind.MOLECULAR,
                           state=LINEAR, inferred=True,
                           source=SOURCE + " (adapter bases proprietary)")
    opposite = graph.product("opposite_strand", "non-retained complementary library strand",
                             MaterialKind.WASTE)
    graph.react("denature", "Denature the indexed double-stranded library",
                StepKind.STRAND_SEPARATION,
                ((ds_library, InputRole.TEMPLATE), (denaturant, InputRole.CONDITION)),
                (strand, opposite))
    graph.discard(opposite, "is not the strand selected for ss-circle preparation")
    splint = graph.input("circle_splint", "proprietary dual-barcode circle splint",
                         MaterialKind.OLIGO, inferred=True,
                         source=SOURCE + " (bases proprietary)")
    ligase = graph.input("circle_ligase", "DNB library circularization ligase",
                         MaterialKind.CATALYST)
    circle = graph.product("ss_circle", "splint-ligated ssDNA library circle",
                           MaterialKind.MOLECULAR, state=CIRCLE, inferred=True,
                           source=SOURCE + " (closure sequences proprietary)")
    graph.react("circularize", "Splint-ligate the library strand into a closed circle",
                StepKind.CIRCULARIZATION,
                ((strand, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.unavailable(
                    "The manual specifies splint-assisted ssDNA closure but does not disclose the aligned adapter and splint bases."))
    exonuclease = graph.input("exonuclease", "linear-DNA exonuclease",
                              MaterialKind.CATALYST)
    selected = graph.product("purified_circle", "exonuclease-purified ssDNA circle",
                             MaterialKind.MOLECULAR, state=CIRCLE, inferred=True,
                             source=SOURCE + " (circle bases proprietary)")
    linear_waste = graph.product("linear_waste", "digested unclosed linear library molecules",
                                 MaterialKind.WASTE)
    graph.react("digest_linear", "Digest remaining linear library molecules",
                StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (selected, linear_waste))
    graph.discard(linear_waste, "is not a covalently closed library circle")
    phi29 = graph.input("phi29", "rolling-circle replication polymerase",
                        MaterialKind.CATALYST)
    dnb = graph.product("dna_nanoball", "DNA nanoball sequencing library",
                        MaterialKind.LIBRARY, state=RCA, inferred=True,
                        source=SOURCE + " (repeat architecture known; adapter bases proprietary)")
    graph.react("rolling_circle", "Copy and compact the circle into a DNA nanoball",
                StepKind.EXTENSION,
                ((selected, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (dnb,))
    graph.sequence(dnb, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="The DNBSEQ manual establishes cPAS entry at the repeated adapter site but keeps the current sequencing-primer bases proprietary."))
    graph.declare_reviewed_scope((
        "indexed_library", "denaturant", "denature", "circle_splint",
        "circle_ligase", "circularize", "exonuclease", "digest_linear", "phi29",
        "rolling_circle", "dna_nanoball"), SOURCE)
    return graph.require_publishable()
