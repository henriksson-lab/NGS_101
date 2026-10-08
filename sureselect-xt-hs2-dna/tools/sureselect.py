"""SureSelect XT HS2 target-capture model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from targeted_ngs import hybrid_capture_scene, sureselect_library

TITLE = "SureSelect XT HS2 DNA target enrichment"
NOTES = "01_sureselect-xt-hs2.html"
SOURCE = 'Agilent <a href="https://www.agilent.com/cs/library/usermanuals/public/G9957-90000.pdf">SureSelect XT HS2 DNA Kits, G9957-90000 B0 (2026)</a>.'
SUMMARY = "A dual-indexed, optionally duplex-MBC-tagged Illumina library is pooled, hybridized to biotinylated RNA baits, recovered on streptavidin beads and amplified after capture."
CAVEAT = "The current manual publishes the molecular-barcode length and index sequences but not every adaptor-oligo base. Kit-specific hidden sequence is represented by its role; no sequence is fabricated."
FINAL_LIBRARY, SEQ_PRIMERS = sureselect_library()
FINAL_CAPTION = "Captured HS2 library in the MBC configuration. Each read begins with its five-base end tag; the two tags identify the original duplex family."
SEQUENCING_INTRO = "The manual specifies standard Illumina paired-end and dual-index primers with no custom sequencing primers."

def sections():
    return [
        ("Prepare and index the library", [
            Row(chunks=[("gDNA → shear → end repair + dA-tail → HS2 adaptor ligation", None, False)]),
            Row(chunks=[("optional MBCs sit immediately beside both insert ends", "umi", False)]),
            Row(chunks=[("pre-capture PCR completes P5 + i5 and i7 + P7", "cbc", False)]),
        ], "Libraries are already sample-indexed before several samples are pooled into one capture."),
        ("Hybridize the target interval to a biotinylated bait",
         hybrid_capture_scene().rows(),
         "Adaptor blockers suppress arm-to-arm hybridization. A design-specific RNA bait selects its complementary insert interval."),
        ("Capture, wash and amplify", [
            Row(chunks=[("bait—library duplex + streptavidin bead → retained", None, False)]),
            Row(chunks=[("unbound library → stringent wash → discarded", None, False)]),
            Row(chunks=[("bead-bound target library → post-capture PCR", None, False)]),
        ], "Affinity selection changes abundance, not adapter orientation or sequencing-primer geometry."),
    ]
