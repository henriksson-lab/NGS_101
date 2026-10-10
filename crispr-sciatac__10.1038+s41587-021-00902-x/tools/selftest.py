#!/usr/bin/env python3
"""Source-transcription checks for CRISPR-sciATAC."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("two combinatorial index parts",P.CELL_PARTS,(("transposition index",8),("PCR-well index",12))); c.report()
