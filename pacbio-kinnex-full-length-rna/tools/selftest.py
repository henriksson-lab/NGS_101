#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import kinnex as K
from checks import Check
c=Check(); c.section("Kinnex source transcription")
c("parallel orientation-specific PCR reactions", K.KINNEX_PCR_REACTIONS, 8)
c("array-formation incubation", K.ARRAY_FORMATION, (45, 60))
c("post-array repair incubation", K.DNA_REPAIR, (45, 30))
c.report()
