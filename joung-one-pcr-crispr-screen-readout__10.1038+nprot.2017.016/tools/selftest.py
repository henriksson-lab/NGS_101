#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
from checks import Check
c = Check()
c("one-PCR library is 260 nt for a non-G-starting guide", len(M.one_pcr_library()), 260)
c.report()
