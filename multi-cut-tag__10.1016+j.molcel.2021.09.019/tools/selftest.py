#!/usr/bin/env python3
"""Source-transcription checks for Multi-CUT&Tag."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("target identity is explicit",[s.feature.role for s in P.FINAL_LIBRARY if s.feature],["sample_index","cell_barcode","feature_barcode","sample_index"]); c.report()
