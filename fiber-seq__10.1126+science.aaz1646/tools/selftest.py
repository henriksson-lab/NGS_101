#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import fiber_seq as M
from checks import Check
c=Check(); c("SMRTbell has no free-end topology", len(M.SMRTBELL.circular_sequence)>len(M.SMRTBELL.insert)*2, True); c.report()
