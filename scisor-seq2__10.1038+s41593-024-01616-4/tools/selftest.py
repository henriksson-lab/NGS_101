#!/usr/bin/env python3
"""Source-transcription checks for ScISOr-seq2."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("Partial Read 1 handle",P.INSERT.segments[0].top,"CTACACGACGCTCTTCCGATCT"); c.report()
