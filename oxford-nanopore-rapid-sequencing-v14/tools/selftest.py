#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import rapid_v14 as M
from checks import Check
c=Check(); c.section("SQK-RAD114 source transcription")
c("fragmentation programme", M.FRAGMENTATION_PROGRAM, ((30, 2), (80, 2)))
c("Rapid Adapter attachment minutes", M.RAPID_ATTACHMENT_MIN, 5)
c.report()
