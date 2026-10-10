#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check
c=Check(); c("MEDS-A transcription", chem.MEDS_A, "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"); c.report()
