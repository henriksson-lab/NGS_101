#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check
import shotgun as M
c=Check()
c.section("Illumina official adapter-table transcription")
for side,printed in (("i5","AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC"),("i7","CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG")):
 parts=list(M.index_primer(side))
 actual=parts[0].top+"["+side+"]"+parts[-1].top
 c(side+" primer arms",actual,printed)
c.report()
