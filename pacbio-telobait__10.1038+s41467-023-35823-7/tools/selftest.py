#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telobait as M
from checks import Check,run_common
c=Check(); c.section("Telobait source choices"); c("six terminal phases",M.R_PHASES,("CCCTAA","ACCCTA","AACCCT","TAACCC","CTAACC","CCTAAC")); c("F01 contains EcoRI", "GAATTC" in M.F01,True); c("capture/release workflow",M.CAPTURE.steps(),("bind 3′ biotin-tagged molecules","wash away untagged DNA","release with EcoRI")); run_common(c); c.report()
