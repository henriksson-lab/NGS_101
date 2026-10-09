"""Oxford Nanopore SQK-RAD114 rapid genomic-DNA workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, annotation_rows


def role(name: str, text: str, tag: str | None = None) -> Segment:
    return Segment(name, text, tag, placeholder=True)


TITLE = "Oxford Nanopore Rapid Sequencing V14 — SQK-RAD114"
NOTES = "01_rapid-sequencing-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/rapid-sequencing-sqk-rad114">'
          'Rapid Sequencing Kit V14 — gDNA, RSE_9177_v114 Rev Q</a>.')
SUMMARY = ("A transposase fragments high-molecular-weight genomic DNA while installing "
           "rapid-attachment tags; the motor-loaded Rapid Adapter then attaches without "
           "ligase and presents one tagged strand to a nanopore.")
CAVEAT = ("The tag and Rapid Adapter bases, attachment chemistry and exact motor-bearing "
          "strand are proprietary. Bracketed regions state only the roles disclosed by "
          "Oxford Nanopore.")
TAGGED = Construct([role("left transposase tag", "[transposase tag]", "me"),
                    role("genomic DNA", "[genomic DNA fragment]"),
                    role("right transposase tag", "[transposase tag]", "me")],
                   name="transposase-tagged genomic DNA")
FINAL_LIBRARY = Construct([
    role("motor-loaded Rapid Adapter", "[motor + Rapid Adapter]", "r2"),
    role("attachment boundary", "[attachment]", "me"),
    role("transposase tag", "[transposase tag]", "me"),
    role("genomic DNA", "[genomic DNA fragment]"),
], name="SQK-RAD114 sequenceable strand")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("No synthesis primer binds. The adapter motor docks at a nanopore "
                     "and meters the attached DNA strand through the pore, adapter-bearing "
                     "end first.")
SEQUENCING_UNAVAILABLE = ("SQK-RAD114 uses motor-controlled nanopore translocation, not "
                          "primer extension; current adapter bases are undisclosed.")
FINAL_CAPTION = ("One sequenceable product. Rapid attachment is non-ligase chemistry, so "
                 "the boundary is labelled attachment rather than covalent ligation.")
FRAGMENTATION_PROGRAM = ((30, 2), (80, 2))
RAPID_ATTACHMENT_MIN = 5


def sections():
    sc = Scene.duplex(list(TAGGED), label="tagged fragment",
                      unpaired=("left transposase tag'", "genomic DNA'",
                                "right transposase tag'"))
    return [
        ("Fragment and tag genomic DNA", sc.rows(),
         "Fragmentation Mix cuts the input and transfers transposase tags in the same "
         "30 °C reaction; 80 °C ends the two-minute tagmentation."),
        ("Attach the Rapid Adapter", [
            Row(chunks=[("[motor + Rapid Adapter]", "r2", False),
                        (" ** attachment ** ", None, False),
                        ("[transposase tag]", "me", False),
                        ("[genomic DNA fragment]", None, False)]),
        ], "Diluted RA attaches to a transposase-tagged end for five minutes at room "
           "temperature without a ligation enzyme."),
        ("Enter the nanopore", [
            Row(chunks=[("membrane ─ pore | motor ─ adapter ─ tagged DNA  ─────────>",
                        None, False)]),
            Row(chunks=[("                     no sequencing primer", "r1", False)]),
        ], "The motor, rather than a polymerase primer, controls strand motion through "
           "the R10.4.1 nanopore."),
    ]
