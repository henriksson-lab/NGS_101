#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import taps
from checks import Check
c=Check(); c.section("source transcription"); c("TruSeq Index 6", taps.INDEX6, "GCCAAT"); c.report()
