#!/usr/bin/env python3
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import chem
import rna_end_readout as shared
from checks import Check
c=Check()
c("Supplementary Figure S1 RT primer", shared.SS2_RT_SEQUENCE,
  "CAGACGTGTGCTCTTCCGATCNNNNNN")
c("custom sequencing primer is the forward-primer suffix",
  shared.SS2_FORWARD_SEQUENCE.endswith(shared.SS2_CUSTOM.top), True)
c.report()
