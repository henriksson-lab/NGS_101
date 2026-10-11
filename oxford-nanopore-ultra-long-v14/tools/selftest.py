#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ultra_long as M
from checks import Check
c=Check(); c.section("SQK-ULK114 source transcription")
c("FRA dilution volumes", M.FRA_DILUTION_UL, (6, 244))
c("tagmentation programme", M.TAGMENTATION_PROGRAM,
  (("room temperature", 10), (75, 10), ("ice", 10)))
c("Rapid Adapter minutes", M.RAPID_ADAPTER_MIN, 30)
c.report()
