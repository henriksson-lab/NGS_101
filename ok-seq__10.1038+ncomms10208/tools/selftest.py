#!/usr/bin/env python3
"""Source-transcription checks for OK-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("four strand-specific adapter oligos",len(P.oligos()),4); c.report()
