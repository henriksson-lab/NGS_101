#!/usr/bin/env python3
"""Source-transcription check for the CG000731 gel-bead oligo."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v4 as V
got="".join(s.top for s in V.bead_oligo()).replace("B","N").replace("U","N")
want="CTACACGACGCTCTTCCGATCT"+"N"*16+"N"*12+"T"*30+"VN"
if got != want: raise SystemExit("FAIL: CG000731 gel-bead oligo transcription")
raise SystemExit(0)
