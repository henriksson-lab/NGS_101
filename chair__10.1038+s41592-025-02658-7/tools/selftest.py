#!/usr/bin/env python3
"""ChAIR: bridge-linker transcription and the computed linker geometry."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from checks import Check

c = Check()
c.section("ChAIR bridge linker (Methods JATS line 179; SI Addendum)")
c("top oligo bases, /iBIOdT/ read as T", chem.BRIDGE_TOP, "CGTGATATTTCACGACTCT")
c("bottom oligo bases", chem.BRIDGE_BOTTOM, "GAGTCGTGAAATATCACGT")
c("annealed ends: 18 bp, 3'-T on both strands", chem.bridge_linker_ends(), (18, "T", "T"))
chem.protocol_graph()
c.report()
