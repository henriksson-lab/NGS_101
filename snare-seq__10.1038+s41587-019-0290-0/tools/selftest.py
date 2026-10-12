#!/usr/bin/env python3
"""Source-transcription checks for SNARE-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("source splint arm", P.SPLINT[:33], "CTGTCTCTTATACACATCTGACGCTGCCGACGA")
c("source splint poly(A)", P.SPLINT[33:], "A" * 31)
c("RNA and accessibility outputs", len(P.FINAL_LIBRARIES), 2)
c.report()
