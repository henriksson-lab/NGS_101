#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import stevens as M
from checks import Check, run_common
c=Check(); c.section("Stevens protocol choices"); c("enzyme",M.WORKFLOW.enzyme.name,"AluI"); c("blunt",M.WORKFLOW.enzyme.end,"blunt"); c("capture",M.WORKFLOW.capture,True); run_common(c); c.report()
