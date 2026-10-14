#!/usr/bin/env python3
"""DSBCapture: the oligo set exists only as figure artwork.

Lensing et al. draw the modified P5/P7 duplexes, the control modified P5 duplex and the
PCR primers in Fig. 1a and Supplementary Fig. 1c; no base of them occurs as text in the
author manuscript, Supplementary Figures or Supplementary Tables.  ``lib/dsb_capture.py``
therefore assembles the set from the shared ``illumina`` constants and derives each bottom
strand and PCR primer from its own top strand, so restating those bases here would only
duplicate a transcription this repository cannot re-check.  What is checked instead is the
one hand-transcribed string, the two end geometries that the derivations do not imply, and
the evidence signalling that the figure-only index window requires.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem
from dsb_capture import DSB_AD006, DSB_P5, DSB_P5_BOTTOM, DSB_P7, DSB_P7_BOTTOM

c = Check()
c.section("DSBCapture Supplementary Fig. 1c transcription")
# The only DSBCapture base string not derived from lib/illumina or another strand.
c("example index AD006 as drawn", DSB_AD006, "GCCAAT")

c.section("Drawn ligation-end geometry")
c("the modified P5 duplex presents one unpaired 3' T",
  (DSB_P5[-1], len(DSB_P5) - len(DSB_P5_BOTTOM)), ("T", 1))
c("the modified P7 duplex presents its 3' T on the bottom strand",
  (DSB_P7_BOTTOM[-1], len(DSB_P7_BOTTOM) - len(DSB_P7)), ("T", 1))

c.section("Figure-only evidence is marked on the page")
inferred = tuple(s.name for s in chem.FINAL_LIBRARY if s.inferred)
c("only the figure-measured index window is shaded inferred", inferred,
  ("6-base i7 reverse complement",))
c("the preamble caveat names the artwork and the index window",
  ("Supplementary Fig. 1c" in chem.CAVEAT, "index" in chem.CAVEAT), (True, True))
graph = chem.protocol_graph()
c("the index-installing step carries an INFERRED caption",
  graph.steps["library_pcr"].label.startswith("INFERRED — "), True)
c("the index read's handling is left undetermined",
  chem.READ_USES, {"Index 1 (i7)": "undetermined"})

sys.path.insert(0, str(HERE.parents[1] / "tools"))
from source_oligo_check import candidates  # noqa: E402
c("every published oligo is a graph state excused by name, with a reason",
  sorted({x.material for x in candidates(graph)}), sorted(chem.SOURCE_OLIGO_EXCEPTIONS))
c.report()
