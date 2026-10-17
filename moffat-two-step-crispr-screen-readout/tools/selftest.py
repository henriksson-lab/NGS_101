#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import build_page as B
import pooled_crispr as M
from chemdraw import revcomp
from checks import Check
c = Check()
c("PCR1 carries no P5", not B.pcr1_product().top().startswith("AATGATAC"))
# The shared library model takes its i5/i7 from lib/illumina; Table 1 prints full Ultramers.
lib = M.two_step_library().top()
c("PCR2 library begins with the full Table 1 D501-F Ultramer", lib.startswith(B.D501_F), True)
c("PCR2 library ends with the complement of the full Table 1 D701-R Ultramer",
  lib.endswith(revcomp(B.D701_R)), True)
c.report()
