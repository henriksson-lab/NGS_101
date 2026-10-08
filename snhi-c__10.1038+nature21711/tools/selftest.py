#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import snhic as M
from checks import Check, run_common
c=Check(); c.section("Flyamer protocol choices"); c("enzyme",M.WORKFLOW.enzyme.name,"DpnII"); c("no biotin",M.WORKFLOW.marking,"none"); c("amplification",M.WORKFLOW.amplification,"MDA"); run_common(c); c.report()
