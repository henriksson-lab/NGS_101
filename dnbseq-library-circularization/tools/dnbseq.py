"""MGI/DNBSEQ ss-circle and DNA-nanoball preparation model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Row, circle_rows, feature
from circular import CircularizationGeometry, circularize_ssdna, rolling_circle
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind, Terminal)

TITLE = "DNBSEQ library circularization and DNA-nanoball preparation"
NOTES = "01_dnbseq.html"
SOURCE = ('Commercial protocol: MGI <a href="https://global-mgitech.com/wp-content/uploads/Documents/SOP-013-B02-185-7.0-MGIEasy-Dual-Barcode-Circularization-Kit-User-Manual-2.pdf">MGIEasy Dual Barcode Circularization Kit User Manual 7.0</a> (part SOP-013-B02-185, kit version V1.0, Sep. 2024) and '
          '<a href="https://en.mgi-tech.com/Uploads/Temp/file/20240530/66583952cae5b.pdf">DNBSEQ Dual Barcode Adapter &amp; Barcode Sequences 3.0</a>.')
SUMMARY = "An indexed double-stranded library is denatured, one strand is splint-ligated into a single-stranded circle, remaining linear DNA is digested, and rolling-circle replication compacts tandem copies into a DNA nanoball for patterned-array sequencing."
CAVEAT = "MGI publishes adapter trimming arms and index choices for its UDB/UMI library-prep family, but the circularization kit accepts any MGI dual-barcode library (Duplex UMI, UDB or &ldquo;MGI Dual-barcode standard product&rdquo;; the manual states it is <i>not</i> suitable for single-barcode libraries). The checked manual supplies reaction volumes, temperatures and QC but names its reagents only as Dual Barcode Splint Buffer, DNA Rapid Ligase, Digestion Buffer and Digestion Enzyme: no splint bases, no enzyme identities and no adapter junction. Its own endpoint is the purified ssCir library; DNA-nanoball formation is drawn here from the platform architecture, not from this manual, and a separate sequencing-primer reference remains unread. The circle topology below therefore uses placeholders."
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
         "Upstream fragmentation, repair, adapter ligation and PCR produce the platform library; 1 pmol of that dual-barcode PCR product (minimum 0.5 pmol) is taken to 48 µL in TE and heat-denatured at 95 °C for 3 min, then chilled on ice."),
        ("Splint-ligate the strand into a closed circle",
         circular.rows()+circle_rows(CIRCLE.linear, closure_label="splint ligation"),
         "12.1 µL of Dual Barcode Splint Buffer plus DNA Rapid Ligase is added and incubated at 37 °C for 30 min; a complementary splint aligns the two adapter ends so the closure joins the final adapter segment back to the first. A 4 µL Digestion Buffer/Digestion Enzyme reaction (37 °C, 30 min, stopped with 7.5 µL Digestion Stop Buffer) then removes molecules that did not close; the manual does not name the enzyme."),
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
    denaturant = graph.input("denaturant", "95 °C, 3 min heat denaturation of 1 pmol PCR product in 48 µL TE, then ice",
                             MaterialKind.CONDITION)
    strand = graph.product("library_strand", LINEAR.name, MaterialKind.MOLECULAR,
                           state=LINEAR, inferred=True,
                           source=SOURCE + " (input-kit-specific adapter structure unresolved)")
    opposite = graph.product("opposite_strand", "complementary library strand",
                             MaterialKind.MOLECULAR)
    graph.react("denature", "Denature the indexed double-stranded library",
                StepKind.STRAND_SEPARATION,
                ((ds_library, InputRole.TEMPLATE), (denaturant, InputRole.CONDITION)),
                (strand, opposite))
    graph.mark_terminal(Terminal(
        opposite.id, Disposition.UNRESOLVED,
        "The manual denatures the dual-barcode PCR product and adds a single splint mixture; it does not state which strand is circularized or what becomes of its complement."))
    splint = graph.input("circle_splint", "Dual Barcode Splint Buffer, 11.6 µL (splint oligo bases undisclosed)",
                         MaterialKind.OLIGO, inferred=True,
                         source=SOURCE + " (splint bases absent from checked manual)")
    ligase = graph.input("circle_ligase", "DNA Rapid Ligase, 0.5 µL; 37 °C for 30 min",
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
    exonuclease = graph.input("digestion_enzyme",
                              "kit Digestion Enzyme (2.6 µL) in Digestion Buffer (1.4 µL); 37 °C for 30 min, stopped with Digestion Stop Buffer",
                              MaterialKind.CATALYST, inferred=True,
                              source=SOURCE + " (enzyme identity and activity not disclosed)")
    selected = graph.product("purified_circle", "exonuclease-purified ssDNA circle",
                             MaterialKind.MOLECULAR, state=CIRCLE, inferred=True,
                             source=SOURCE + " (input-kit-specific circle bases unresolved)")
    linear_waste = graph.product("linear_waste", "digested non-circular library molecules",
                                 MaterialKind.WASTE)
    graph.react("digest_linear", "Digest the molecules that did not circularize",
                StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (selected, linear_waste))
    graph.discard(linear_waste, "is not a covalently closed library circle; the kit's bead cleanup then removes the digestion products")
    beads = graph.input("clean_beads",
                        "MGIEasy DNA Clean Beads, 170 µL, two 80% ethanol washes, eluted in 22 µL TE (20 µL recovered)",
                        MaterialKind.SELECTION)
    pure_circle = graph.product("clean_circle", "purified ssCir library",
                                MaterialKind.MOLECULAR, state=CIRCLE, inferred=True,
                                source=SOURCE + " (input-kit-specific circle bases unresolved)")
    bead_waste = graph.product("cleanup_waste", "digestion reagents and short products in the discarded supernatant and washes",
                               MaterialKind.WASTE)
    graph.react("clean_up", "Bead-clean the digestion product", StepKind.CLEANUP,
                ((selected, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (pure_circle, bead_waste))
    graph.discard(bead_waste, "is removed with the supernatant and ethanol washes")
    phi29 = graph.input("phi29", "rolling-circle replication polymerase (DNB-preparation reagent; not part of the circularization kit)",
                        MaterialKind.CATALYST)
    dnb = graph.product("dna_nanoball", "DNA nanoball sequencing library",
                        MaterialKind.LIBRARY, state=RCA, inferred=True,
                             source=SOURCE + " (the circularization manual stops at the purified ssCir library; rolling-circle replication and nanoball loading come from the DNBSEQ platform architecture, not from this manual)")
    graph.react("rolling_circle", "Copy and compact the circle into a DNA nanoball",
                StepKind.EXTENSION,
                ((pure_circle, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)), (dnb,))
    graph.sequence(dnb, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="The checked DNBSEQ manuals establish cPAS entry at the repeated adapter site but do not give the sequencing-primer bases; the separate primer reference guide remains unreviewed."))
    graph.annotate(Highlight(
        "clean_circle", "The checked kit manual ends here",
        "Manual 7.0 (sections 3.1–3.4) covers denaturation, splint ligation, digestion, bead cleanup and Qubit ssDNA QC of the purified ssCir library, with a specified yield of at least 80 fmol. It contains no DNB, rolling-circle or sequencing-primer text; the steps drawn after this point describe the DNBSEQ platform rather than this kit, and their reagents and bases remain unverified here."))
    graph.annotate(Highlight(
        "indexed_library", "Only dual-barcode libraries are accepted",
        "Sections 1.1–1.2 restrict the input to PCR products carrying MGI Duplex UMI or UDB adapters, or an &ldquo;MGI Dual-barcode standard product&rdquo;, and state that the kit is not suitable for circularizing single-barcode libraries. The recommended input is 1 pmol, with 0.5 pmol as the stated minimum."))
    graph.declare_reviewed_scope((
        "indexed_library", "denaturant", "denature", "opposite_strand", "circle_splint",
        "circle_ligase", "circularize", "digestion_enzyme", "digest_linear",
        "clean_beads", "clean_up", "clean_circle", "phi29",
        "rolling_circle", "dna_nanoball"), SOURCE)
    return graph.require_publishable()
