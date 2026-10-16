#!/usr/bin/env python3
"""DSBCapture published index and ligation geometry checks.

The authors' original Protocol Exchange method prints the oligos and six-base
single-index design. tools/selftest_source_oligos.py compares the graph's exact
oligo states with that source when the source bundle is available.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem
from chem import DSB_AD006, DSB_P5, DSB_P5_BOTTOM, DSB_P7, DSB_P7_BOTTOM

c = Check()
c.section("DSBCapture original detailed method")
c("AD06 index agrees with the printed PCR R example", DSB_AD006, "GCCAAT")

c.section("Drawn ligation-end geometry")
c("the modified P5 duplex presents one unpaired 3' T",
  (DSB_P5[-1], len(DSB_P5) - len(DSB_P5_BOTTOM)), ("T", 1))
c("the modified P7 duplex presents its 3' T on the bottom strand",
  (DSB_P7_BOTTOM[-1], len(DSB_P7_BOTTOM) - len(DSB_P7)), ("T", 1))

c.section("Published evidence is marked on the page")
inferred = tuple(s.name for s in chem.FINAL_LIBRARY if s.inferred)
c("published index window is established", inferred, ())
c("the preamble cites the detailed method", "protex.2016.52" in chem.CAVEAT, True)
graph = chem.protocol_graph()
c("the PCR step states the published index width",
  "six-base" in graph.steps["library_pcr"].label, True)
c("the index read's handling is left undetermined",
  chem.READ_USES, {"Index 1 (i7)": "undetermined"})

c.report()
