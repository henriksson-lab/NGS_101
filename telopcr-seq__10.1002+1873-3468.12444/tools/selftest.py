#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telopcr as M
from checks import Check,run_common
c=Check(); c.section("TeloPCR-seq source choices"); c("S. pombe repeat",M.END.repeat,"GGTTACA"); c("anchor remains explicit unknown sequence",set(M.ANCHOR),{"N"}); run_common(c); c.report()
