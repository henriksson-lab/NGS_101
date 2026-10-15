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
SOURCE = ('Commercial protocol: MGI <a href="https://global-mgitech.com/wp-content/uploads/Documents/SOP-013-B02-185-7.0-MGIEasy-Dual-Barcode-Circularization-Kit-User-Manual-2.pdf">Dual Barcode Circularization Kit Manual 7.0</a> and '
          '<a href="https://en.mgi-tech.com/Uploads/Temp/file/20240530/66583952cae5b.pdf">DNBSEQ Dual Barcode Adapter &amp; Barcode Sequences 3.0</a>.')
SUMMARY = "An indexed double-stranded library is denatured, one strand is splint-ligated into a single-stranded circle, remaining linear DNA is digested, and rolling-circle replication compacts tandem copies into a DNA nanoball for patterned-array sequencing."
CAVEAT = "MGI publishes adapter trimming arms and index choices for its UDB/UMI library-prep family, but the circularization kit accepts libraries from different preparation routes. The checked circularization manual does not give the splint sequence or establish an exact adapter junction for every input kit; a separate sequencing-primer reference remains to be checked. The generic circle topology below therefore uses placeholders."
SEQ_PRIMERS = ()
FINAL_LIBRARY = None
SEQUENCING_ENDING = "A DNBSEQ sequencing primer anneals to the repeated adapter site within each DNA nanoball. Because rolling-circle replication creates many tandem copies of the same ssDNA circle, every repeat presents the same primer site and insert boundary to cPAS sequencing. The checked manuals do not give the sequencing-primer bases, so the binding geometry uses a placeholder."

def oligos():
    return [
        "MGI UDB/UMI prep-family adapter trimming arms (manufacturer reference 3.0; not a complete circularization input):",
        "Forward: 5′-AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA-3′",
        "Reverse: 5′-AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG-3′",
    ]

LINEAR = Construct([seg("adapter A / primer site", "X" * 16, "r1", placeholder=True,
                        inferred=True),
                    seg("insert", "N" * 30, placeholder=True),
                    seg("adapter B / barcode", "Y" * 16, "cbc", placeholder=True,
                        inferred=True, feature=feature("sample_barcode_region", "sample_index", "unknown",
                                                       note="barcode sub-boundary depends on the selected input kit"))], name="denatured DNBSEQ library strand")
CIRCLE = circularize_ssdna(LINEAR, five_prime_phosphate=True)
RCA = rolling_circle(CIRCLE, "X" * 16, copies=3)

def sections():
    upstream=AdapterLigationGeometry.unavailable("The circularization kit accepts several upstream library-prep inputs; the checked UDB/UMI adapter trimming sequences do not establish every ligation-facing end geometry.")
    circular=AdapterLigationGeometry.unavailable("The DNB-preparation manual describes splint-assisted circle closure but does not disclose the splint bases or the aligned adapter-end sequences.")
    return [
        ("Denature the indexed double-stranded library",
         upstream.rows()+[Row(chunks=[("indexed dsDNA → heat/alkali → one adapter-bearing ssDNA strand", None, False)])],
         "Upstream fragmentation, repair, adapter ligation and PCR produce the platform library; DNB preparation begins by making it single-stranded."),
        ("Splint-ligate the strand into a closed circle",
         circular.rows()+circle_rows(CIRCLE.linear, closure_label="splint ligation"),
         "A complementary splint aligns the two adapter ends for ligation. The closure joins the final adapter segment back to the first; exonuclease removes unclosed linear molecules."),
        ("Copy the circle into a compact DNA nanoball",
         RCA.rows("rolling-circle product before compaction")+
         [Row(chunks=[("self-compaction → one DNA nanoball", None, False)])],
         "Rolling-circle replication produces a single concatemer containing tandem complementary copies of the original circle; it is compacted and loaded onto a patterned array."),
        ("Prime repeated adapter sites for cPAS sequencing",
         [Row(chunks=[("…[primer site]—insert—barcode—[primer site]—insert—barcode…", "r1", False)]),
          Row(chunks=[("   sequencing primer →              sequencing primer →", "r2", False)])],
         "The same as-yet-unresolved sequencing-primer site occurs once per RCA repeat. Probe-anchor synthesis interrogates the immobilized nanoball."),
    ]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    ds_library = graph.input("indexed_library", "adapter-bearing indexed dsDNA library",
                             MaterialKind.MOLECULAR, primary=True)
    denaturant = graph.input("denaturant", "heat or alkaline denaturation condition",
                             MaterialKind.CONDITION)
    strand = graph.product("library_strand", LINEAR.name, MaterialKind.MOLECULAR,
                           state=LINEAR, inferred=True,
                           source=SOURCE + " (input-kit-specific adapter structure unresolved)")
    opposite = graph.product("opposite_strand", "non-retained complementary library strand",
                             MaterialKind.WASTE)
    graph.react("denature", "Denature the indexed double-stranded library",
                StepKind.STRAND_SEPARATION,
                ((ds_library, InputRole.TEMPLATE), (denaturant, InputRole.CONDITION)),
                (strand, opposite))
    graph.discard(opposite, "is not the strand selected for ss-circle preparation")
    splint = graph.input("circle_splint", "dual-barcode circle splint (sequence unresolved)",
                         MaterialKind.OLIGO, inferred=True,
                         source=SOURCE + " (splint bases absent from checked manual)")
    ligase = graph.input("circle_ligase", "DNB library circularization ligase",
                         MaterialKind.CATALYST)
    circle = graph.product("ss_circle", "splint-ligated ssDNA library circle",
                           MaterialKind.MOLECULAR, state=CIRCLE, inferred=True,
                           source=SOURCE + " (input-kit-specific closure sequences unresolved)")
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
                             source=SOURCE + " (input-kit-specific circle bases unresolved)")
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
                             source=SOURCE + " (repeat architecture known; input adapter structure unresolved)")
    graph.react("rolling_circle", "Copy and compact the circle into a DNA nanoball",
                StepKind.EXTENSION,
                ((selected, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (dnb,))
    graph.sequence(dnb, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="The checked DNBSEQ manuals establish cPAS entry at the repeated adapter site but do not give the sequencing-primer bases; the separate primer reference guide remains unreviewed."))
    graph.declare_reviewed_scope((
        "indexed_library", "denaturant", "denature", "circle_splint",
        "circle_ligase", "circularize", "exonuclease", "digest_linear", "phi29",
        "rolling_circle", "dna_nanoball"), SOURCE)
    return graph.require_publishable()
