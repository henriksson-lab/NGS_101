#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import trac_looping as P
from checks import Check
c=Check(); c.section("TrAC-looping source transcription")
c("Bio67F length", len(P.BIO67F), 67)
c("67bpR length", len(P.BP67R), 67)
c.report()
