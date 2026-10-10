#!/usr/bin/env python3
"""Source-transcription checks for OAK."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check(); c("droplet and aliquot identity parts",P.CELL_PARTS,(("droplet barcode",16),("aliquot index",8))); c.report()
