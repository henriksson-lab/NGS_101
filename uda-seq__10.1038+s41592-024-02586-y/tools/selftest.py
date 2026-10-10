#!/usr/bin/env python3
"""Source-transcription checks for UDA-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("droplet and well identity parts",P.CELL_PARTS,(("droplet barcode",16),("well barcode",8))); c.report()
