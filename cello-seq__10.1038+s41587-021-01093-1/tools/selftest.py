#!/usr/bin/env python3
"""Source-transcription checks for CELLO-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("patterned UMI length",next(len(s.top) for s in P.FINAL_LIBRARY if s.feature and s.feature.role=="umi"),22); c.report()
