#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from checks import Check

c = Check()
c.section("TIP-seq source transcription")
c("ME-T7-bulk", chem.ME_T7_BULK,
  "GAATTTAATACGACTCACTATAGGGTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("ME-scT7-r5001", chem.ME_SC_R5001,
  "GAATTTAATACGACTCACTATAGGGAGATCCACGCGCTCTAATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("ME reverse", chem.ME_REV, "CTGTCTCTTATACACATCT")
c("second-strand primer", chem.SSS_SCI, "GGGAGATCCACGCGC")
chem.workflow()
c.report()
