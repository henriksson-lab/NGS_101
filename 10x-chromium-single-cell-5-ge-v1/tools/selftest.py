#!/usr/bin/env python3
"""Source-transcription check for the CG000109 gel-bead TSO."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v1 as V
got="".join(s.top for s in V.bead_tso()).replace("B","N").replace("U","N")
want="CTACACGACGCTCTTCCGATCT"+"N"*16+"N"*10+"TTTCTTATATGGG"
if got != want: raise SystemExit("FAIL: CG000109 gel-bead TSO transcription")
raise SystemExit(0)
