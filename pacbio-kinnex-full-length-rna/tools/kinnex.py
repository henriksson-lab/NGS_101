"""PacBio Kinnex full-length RNA ordered-array workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import (Construct, Row, Segment, feature, feature_rows, revcomp,
                      strand_row)
from concatemer import SegmentedArray
from dumbbell import Dumbbell, dumbbell_rows
from rna_special import template_switch_scene


def role(name: str, width: int, tag: str | None = None, **kw) -> Segment:
    return Segment(name, "X" * width, tag, placeholder=True, **kw)


TITLE = "PacBio Kinnex full-length RNA — eight-insert arrays"
NOTES = "01_kinnex-full-length-rna.html"
SOURCE = ('Commercial defining source: PacBio <a href="https://www.pacb.com/wp-content/'
          'uploads/Procedure-checklist-Preparing-Kinnex-libraries-using-the-Kinnex-full-'
          'length-RNA-kit.pdf">Preparing Kinnex libraries using the Kinnex full-length '
          'RNA kit, 103-238-700 Rev09</a>.')
SUMMARY = ("Full-length cDNA is copied in eight parallel orientation-specific PCRs; the "
           "products assemble into ordered eight-insert arrays whose terminal adapters "
           "form nuclease-resistant SMRTbell templates for HiFi sequencing.")
CAVEAT = ("Kinnex primer, segmentation-junction and terminal-adapter bases are proprietary. "
          "Their experimentally specified order and topology are shown as named unknown DNA.")

ISOSEQ_INDEX = feature("isoseq_sample_index", "sample_index", "whitelist",
                       whitelist="PacBio Iso-Seq primer barcodes bc01–12")
INSERTS = tuple(Construct([
                    Segment("Iso-Seq primer barcode", "[Iso-Seq barcode]", "cbc",
                            placeholder=True, inferred=True, feature=ISOSEQ_INDEX),
                    role("full-length cDNA", 18, "r1")],
                          name=f"orientation-specific Kinnex product {i}")
                for i in range(1, 9))
JUNCTIONS = tuple(role(f"segmentation junction {i}|{i + 1}", 8, "r3")
                  for i in range(1, 8))
ARRAY = SegmentedArray(INSERTS, JUNCTIONS,
                       role("left terminal adapter", 12, "r2"),
                       role("right terminal adapter", 12, "r2"),
                       name="eight-insert Kinnex array")
SMRTBELL = Dumbbell(ARRAY.core(), "N" * 18, "N" * 18,
                    name="Kinnex array SMRTbell", insert_overhang="A",
                    adapter_overhang="T")
KINNEX_PCR_REACTIONS = 8
ARRAY_FORMATION = (45, 60)
DNA_REPAIR = (45, 30)
NUCLEASE = (37, 15)


def polymerase_cycle() -> Construct:
    core = ARRAY.core()
    return Construct([
        role("eight cDNA inserts and seven segmentation junctions, forward", len(core), "r1"),
        role("right terminal adapter / primer site", 18, "r2"),
        role("eight cDNA inserts and seven segmentation junctions, reverse", len(core), "r1"),
        role("left terminal adapter / primer site", 18, "r2"),
    ], name="one Kinnex SMRTbell polymerase circuit")


def sections():
    return [
        ("Make full-length cDNA", template_switch_scene(inferred=True).rows(),
         "INFERRED — the current kit discloses the poly(A)-primed, template-switch workflow "
         "but not the oligo bases. Iso-Seq barcode primers are added during cDNA PCR."),
        ("Install eight orientation-specific segmentation ends", [
            Row(chunks=[("one cDNA pool → PCR A | B | C | D | E | F | G | HQ", "r3", False)]),
            strand_row(INSERTS[0]),
            *feature_rows(INSERTS[0], prefix_width=4),
        ], "Eight parallel Kinnex PCRs give each copy the end identities needed for ordered "
           "array formation; the Iso-Seq primer barcode retains sample identity, and equal "
           "volumes of all eight products are pooled."),
        ("Assemble the eight-insert array", ARRAY.rows(),
         "Kinnex enzyme and ligase assemble the eight products between barcoded terminal "
         "adapters. Every ** is a covalent insert–segmentation-junction boundary."),
        ("Select complete SMRTbells with nuclease", dumbbell_rows(SMRTBELL),
         "The terminal adapters close the ordered array. Nuclease removes molecules with "
         "unprotected ends; a complete closed template survives."),
    ]
