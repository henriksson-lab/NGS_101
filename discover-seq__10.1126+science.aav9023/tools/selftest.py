#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from editing_outcomes import DISCOVER_AR001, discover_library

c = Check()
c.section("DISCOVER-Seq source transcription")
c("TruSeq ChIP Table 5 AR001 index", DISCOVER_AR001, "ATCACG")
lib, _ = discover_library()
c("historical library carries a six-base sample index",
  len(next(s for s in lib if s.name == "i7 reverse complement")), 6)
c.report()
