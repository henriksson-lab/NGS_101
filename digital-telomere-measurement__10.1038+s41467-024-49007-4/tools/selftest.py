#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
from checks import Check,run_common
c=Check(); c.section("DTM source boundary"); c("unverified capture core stays unknown",set(M.CAPTURE.core),{"N"}); c("all telomere phases represented",len(M.CAPTURE.arms),6); run_common(c); c.report()
