#!/usr/bin/env python3
"""Source-transcription checks for SUM-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("sample and droplet identity parts",P.CELL_PARTS,(("sample index",8),("droplet barcode",16))); c.report()
