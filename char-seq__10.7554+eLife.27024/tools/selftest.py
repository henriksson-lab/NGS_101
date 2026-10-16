#!/usr/bin/env python3
"""Source-transcription checks for ChAR-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("bridge top transcription", P.BRIDGE_TOP,
  "AAAAAAACCGGCGTCCAAGGATCTAG")
c("bridge bottom transcription", P.BRIDGE_BOTTOM,
  "GATCCTAGATCCTTGGACGCCGGT")
c("six-base unpaired RNA-ligation tail", P.TOP_TAIL, "AAAAAA")
c("2018 bridge has no UMI",
  [s.feature.role for s in P.FINAL_LIBRARY if s.feature].count("umi"), 0)
c("defining run is single-end", P.READ_LENGTHS, {"Read 1": 152})
c.report()
