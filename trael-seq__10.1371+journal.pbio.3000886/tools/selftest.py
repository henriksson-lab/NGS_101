#!/usr/bin/env python3
"""Source-transcription checks for TrAEL-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("two source adapters",len(P.oligos()),2); c.report()
