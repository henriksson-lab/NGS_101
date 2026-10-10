#!/usr/bin/env python3
"""Source-transcription checks for CapTrap-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("native-end selection is represented",[s.name for s in P.FINAL_LIBRARY][1:],["5-prime linker","cap-selected full-length cDNA","copied poly(A) tail","3-prime linker"]); c.report()
