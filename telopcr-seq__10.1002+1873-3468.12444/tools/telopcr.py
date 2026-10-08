"""TeloPCR-seq molecular workflow (Bennett et al. 2016)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from dumbbell import Dumbbell, dumbbell_rows
from telomere import TelomereEnd

END = TelomereEnd(repeat="GGTTACA", duplex_copies=3, overhang_copies=2)
ANCHOR = "N" * 24  # structure published; exact ordered bases are in the inaccessible supplement
AMPLICON = Construct([Segment("subtelomere primer flank", "X" * 18, placeholder=True),
                      Segment("telomere", END.repeat * 4, "r1"),
                      Segment("terminal anchor", ANCHOR, "r3", placeholder=True)], name="TeloPCR amplicon")
SMRTBELL = Dumbbell(AMPLICON, "N" * 16, "N" * 16, name="PacBio circular library")

def panels():
    return [("Attach a single-stranded anchor to the native chromosome terminus", END.scene(),
             "CircLigase joins the anchor to the native telomere end; the exact supplement oligo is not fabricated here."),
            ("Amplify from a subtelomere-specific primer to the anchor primer",
             [Row(chunks=[("subtelomere primer  --->  [subtelomere][telomere repeats][terminal anchor]  <--- anchor primer", None, False)])],
             "Each chromosome end is amplified with its subtelomeric primer and the shared anchor primer."),
            ("Convert the amplicon into a circular PacBio template", dumbbell_rows(SMRTBELL),
             "Hairpin ligation closes both amplicon ends; ** marks the sealed insert–hairpin junctions.")]
