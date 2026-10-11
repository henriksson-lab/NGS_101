#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from bulk_atac import ATAC_AD1_NOMX, ATAC_AD2_1
from checks import Check

c = Check()
c.section("original ATAC amplification-primer transcription")
c("Ad1_noMX", ATAC_AD1_NOMX,
  "AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG")
c("Ad2.1_TAAGGCGA", ATAC_AD2_1,
  "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT")
c.report()
