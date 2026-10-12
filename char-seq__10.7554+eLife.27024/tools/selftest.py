#!/usr/bin/env python3
"""Source-transcription checks for ChAR-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("bridge top transcription", P.BRIDGE_TOP,
  "AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG")
c("bridge bottom transcription", P.BRIDGE_BOTTOM,
  "GATCTGCGACTTAATTAAAGATCCTTGGACGCCGGTT")
c("bridge strands share the exact 31-bp core", P.TOP_CORE,
  "CCGGCGTCCAAGGATCTTTAATTAAGTCGCA")
c("molecular identifier is retained",
  [s.feature.role for s in P.FINAL_LIBRARY if s.feature].count("umi"), 1)
c("defining run is single-end", P.READ_LENGTHS, {"Read 1": 152})
c.report()
