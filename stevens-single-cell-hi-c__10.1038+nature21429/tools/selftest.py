#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import stevens as M
from checks import Check, run_common
c=Check(); c.section("Stevens source transcription"); c("contact enzyme",M.WORKFLOW.enzyme.name,"MboI"); c("cohesive contact end",M.WORKFLOW.enzyme.end,"5-prime"); c("capture",M.WORKFLOW.capture,True); run_common(c); c.report()
