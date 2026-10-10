#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import methyl_hic as M
from checks import Check
c=Check(); c.section("Methyl-HiC source facts")
c("restriction enzyme", M.WORKFLOW.enzyme.name, "DpnII")
c("captured biotin junction", (M.WORKFLOW.marking, M.WORKFLOW.capture), ("biotin-fill", True))
c.report()
