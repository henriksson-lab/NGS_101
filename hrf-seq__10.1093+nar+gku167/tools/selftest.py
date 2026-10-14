#!/usr/bin/env python3
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import rna_end_readout as shared
from checks import Check
c=Check()
c("published random RT primer", shared.HRF_RT_RANDOM,
  "AGACGTGTGCTCTTCCGATCTNNNNNNNNS")
c("four published sample-index oligos", shared.HRF_INDEXES,
  ("GGAACT", "GGACGG", "CGTACG", "GCTACC"))
c.report()
