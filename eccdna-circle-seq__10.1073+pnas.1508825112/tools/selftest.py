#!/usr/bin/env python3
"""Source transcription of the six-base index reagent inventory."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import eccdna_circle as M
from checks import Check,Source,have
from storage import data_root
c=Check()
source=Source(data_root()/"sources"/H.parent.name/"pnas.1508825112_PMC4475933.html.txt", "defining PNAS article Methods Sequencing")
if have(c,source,label="Apollo index-oligo length source"):
    match=re.search(r"barcode index labels \((\d+)-base oligos\)",source.path.read_text())
    inventory=M.apollo_reagent_inventory().strands["barcodes"].segs
    c("source barcode oligo length; no placement inferred",sum(len(s) for s in inventory),int(match.group(1)) if match else None)
c.report()
