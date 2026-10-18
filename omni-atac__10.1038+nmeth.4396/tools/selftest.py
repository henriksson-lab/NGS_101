#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE),str(HERE.parents[1] / "lib")]
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from chem import ATAC_AD1_NOMX, ATAC_AD2_1, protocol_graph
from checks import Check

c = Check()
c.section("original ATAC amplification-primer transcription")
c("Ad1_noMX", ATAC_AD1_NOMX,
  "AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG")
c("Ad2.1_TAAGGCGA", ATAC_AD2_1,
  "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT")
protocol_graph()
c.report()
