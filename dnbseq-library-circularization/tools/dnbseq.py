"""MGI/DNBSEQ ss-circle and DNA-nanoball preparation model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, circle_rows
from circular import circularize_ssdna, rolling_circle

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
                        inferred=True)], name="denatured DNBSEQ library strand")
CIRCLE = circularize_ssdna(LINEAR, five_prime_phosphate=True)
RCA = rolling_circle(CIRCLE, "X" * 16, copies=3)

def sections():
    return [
        ("Denature the indexed double-stranded library",
         [Row(chunks=[("indexed dsDNA → heat/alkali → one adapter-bearing ssDNA strand", None, False)])],
         "Upstream fragmentation, repair, adapter ligation and PCR produce the platform library; DNB preparation begins by making it single-stranded."),
        ("Splint-ligate the strand into a closed circle",
         circle_rows(CIRCLE.linear, closure_label="splint ligation"),
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
