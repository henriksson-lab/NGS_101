#!/usr/bin/env python3
"""Source-transcription check for the CG000108 gel-bead oligo."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v2 as V
from checks import Check
check=Check()
got="".join(s.top for s in V.bead_oligo())
want="CTACACGACGCTCTTCCGATCT"+"N"*16+"N"*10+"T"*30+"VN"
check("CG000108 gel-bead oligo transcription",
      got.replace("B","N").replace("U","N"), want)
check.report()
