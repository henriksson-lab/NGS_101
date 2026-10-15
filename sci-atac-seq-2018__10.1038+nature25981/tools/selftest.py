#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE),
                str(HERE.parents[1] / "sci-atac-seq__10.1126+science.aab1601" / "tools"),
                str(HERE.parents[1] / "lib")]

import chem
import illumina as il
import nextera as nx
from checks import Check

check = Check()
# Cusanovich et al. 2018 Supplementary Table S12, P5_1_PCR_Primer.
check("P5_1_PCR_Primer", il.P5 + "CTCCATCGAG" + nx.S5,
      chem.SOURCE_P5_1)
check("the published ten-base PCR indices are established in the construct",
      [(s.name, s.inferred) for s in chem.FINAL_LIBRARY
       if s.name in ("PCR i5", "PCR i7'")],
      [("PCR i5", False), ("PCR i7'", False)])
check.report()

