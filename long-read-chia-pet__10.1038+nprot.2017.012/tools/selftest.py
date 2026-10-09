#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import long_read_chia_pet as M
from checks import Check
c=Check(); c.section("long-read ChIA-PET source transcription")
c("Bridge linker-F", M.BRIDGE_F, "CGCGATATCTTATCTGACT")
c("Bridge linker-R", M.BRIDGE_R, "GTCAGATAAGATATCGCGT")
c("Tn5 incubation", M.TAGMENTATION, (55, 5))
c("stated PCR-cycle ceiling", M.MAX_PCR_CYCLES, 13)
c.report()
