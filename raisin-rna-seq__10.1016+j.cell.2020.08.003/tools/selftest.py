#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H.parents[1]/"lib")]
from checks import Check
from rare_cell_and_affinity import RAISIN_CST
c=Check(); c.section("RAISIN source transcription")
c("published CST buffer", RAISIN_CST,
  (("CHAPS", "0.49% w/v"), ("NaCl", "146 mM"), ("CaCl2", "1 mM"),
   ("MgCl2", "21 mM"), ("Tris pH 8.0", "10 mM")))
c.report()
