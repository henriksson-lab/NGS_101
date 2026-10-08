#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
from checks import Check
c = Check()
c("PCR1 carries no P5", not M.two_step_intermediate().top().startswith("AATGATAC"))
c("PCR2 carries P5 and P7", M.two_step_library().top().startswith("AATGATAC") and M.two_step_library().top().endswith("ATCTCGTATGCCGTCTTCTGCTTG"))
c.report()
