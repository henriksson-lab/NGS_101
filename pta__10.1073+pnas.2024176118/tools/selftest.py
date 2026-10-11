#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
from checks import Check

c = Check()
c.section("PTA source transcription")
c("undisclosed primer stays undisclosed", P.pta.amplification_product().segments[0].top,
  "[random primer]")
c("terminated product records alpha-thio-ddN chemistry",
  "alpha-thio-ddN" in P.pta.amplification_product().segments[-1].note, True)
c.report()
