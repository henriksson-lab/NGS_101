#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import nagano as M
from checks import Check, run_common
c=Check(); c.section("Nagano protocol choices"); c("enzyme",M.WORKFLOW.enzyme.name,"BglII"); c("isolate after ligation","after bulk proximity ligation" in M.WORKFLOW.isolate,True); c("capture",M.WORKFLOW.capture,True); run_common(c); c.report()
