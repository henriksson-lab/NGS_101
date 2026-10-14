#!/usr/bin/env python3
"""Source-transcription checks for iMARGI."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("two linker random bases are not called a UMI",
  "umi" in [s.feature.role for s in P.FINAL_LIBRARY if s.feature], False)
c("author linker retains phosphorothioate and internal biotin",
  "TCC*GA/iBiodT/CT" in P.oligos()[0], True)
c.report()
