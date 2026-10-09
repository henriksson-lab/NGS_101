"""Illumina Nextera Mate Pair library construction."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Row
from chemdraw import Scene
from mate_pair import (DUPLICATE_JUNCTION, circularization_rows, recovered_truseq_library,
                       recovery_scene, tagmentation_scene)

TITLE = "Nextera Mate Pair — circular junction recovery"
NOTES = "01_nextera-mate-pair.html"
SOURCE = ('Commercial defining source: Illumina, <a href="https://support.illumina.com/'
          'content/dam/illumina-support/documents/documentation/chemistry_documentation/'
          'samplepreps_nextera/nexteramatepair/nextera-mate-pair-reference-guide-15035209-02.pdf">'
          '<i>Nextera Mate Pair Library Prep Reference Guide</i>, document 15035209 v02</a>; '
          '<a href="https://support.illumina.com/content/dam/illumina-marketing/documents/'
          'products/technotes/technote_nextera_matepair_data_processing.pdf">junction-sequence technical note</a>.')
SUMMARY = ("Long genomic fragments are Tn5-tagged, gap-filled and circularized. After shearing, "
           "streptavidin selects short fragments spanning the old circle-closing bond; TruSeq "
           "adapters turn those distant genomic ends into an outward-facing read pair.")
CAVEAT = "The drawn duplicated junction is the common product; the vendor also documents single-junction products."

FINAL_LIBRARY, SEQ_PRIMERS = recovered_truseq_library()
FINAL_CAPTION = ("Single-index TruSeq library containing a selected duplicated mate-pair junction. "
                 "The two genomic pieces flanking it were distant ends of one long input fragment.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, six-base Index 1 and Read 2 primers bind the final "
                    "library. The genomic mates map outward because circularization inverted their relationship.")


def adapter_ligation_scene() -> Scene:
    sc = Scene.duplex(list(FINAL_LIBRARY), label="adapter-ligated captured fragment")
    sc.junction("top", "Read 1 arm", "left distant end", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc


def sections():
    return [
        ("Tagment long genomic DNA", [*tagmentation_scene().rows(),
            Row(chunks=[("strand-displacement fill repairs the nine-base Tn5 gaps", None, False)]),
        ], "Mate Pair Tagment Enzyme fragments and tags long DNA; ** marks the transferred tag boundaries."),
        ("Circularize long fragments", circularization_rows(),
         "Blunt intramolecular ligation closes the two distant ends; exonuclease removes molecules that remain linear."),
        ("Shear and recover the old circle junction", recovery_scene().rows(),
         "Covaris shearing yields short pieces. Streptavidin beads retain the biotin-marked junction-spanning subset."),
        ("End-repair, dA-tail and ligate TruSeq index adapter", adapter_ligation_scene().rows(),
         "The final library uses the TruSeq DNA LT single-index workflow, not Nextera sequencing primers."),
    ]
