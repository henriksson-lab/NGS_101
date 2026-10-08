#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
from checks import Check,run_common
c=Check(); c.section("Telomere Profiling source boundary"); c("unverified TeloTag core stays unknown",set(M.TELOTAG.core),{"N"}); c("capture has release step",len(M.CAPTURE.steps()),3); run_common(c); c.report()
