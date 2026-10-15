#!/usr/bin/env python3
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import chem
import chem
from checks import Check
c=Check()
c("Supplementary Figure S1 RT primer", chem.SS2_RT_SEQUENCE,
  "CAGACGTGTGCTCTTCCGATCNNNNNN")
c("custom sequencing primer is the forward-primer suffix",
  chem.SS2_FORWARD_SEQUENCE.endswith(chem.SS2_CUSTOM.top), True)
c.report()
