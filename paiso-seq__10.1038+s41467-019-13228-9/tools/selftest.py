#!/usr/bin/env python3
"""Source-transcription checks for PAIso-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("TSO-derived handle",P.INSERT.segments[0].top,"AAGCAGTGGTATCAACGCAGAGT"); c.report()
