#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import smrtbell2 as M
from checks import Check
c=Check(); c.section("SMRTbell Express 2.0 source transcription")
c("overhang removal", M.OVERHANG_REMOVAL, (37, 15))
c("damage repair", M.DAMAGE_REPAIR, (37, 30))
c("end prep", M.END_PREP_PROGRAM, ((20, 10), (65, 30)))
c("adapter ligation", M.LIGATION, (20, 60))
c.report()
