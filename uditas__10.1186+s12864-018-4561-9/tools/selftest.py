#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import UDITAS_TOP_N501, UDITAS_BOTTOM
c = Check(); c.section("UDiTaS supplementary oligo transcription")
c("N501 top and bottom", (UDITAS_TOP_N501, UDITAS_BOTTOM),
  ("AATGATACGGCGACCACCGAGATCTACACTAGATCGCNNNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "CTGTCTCTTATACAC"))
c.report()
