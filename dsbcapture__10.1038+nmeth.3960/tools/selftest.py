#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from dsb_capture import DSB_P5, DSB_P7, DSB_PCR_F
c = Check()
c.section("DSBCapture Figure 1 transcription")
c("modified P5/P7 and PCR F", (DSB_P5, DSB_P7, DSB_PCR_F),
  ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
   "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
   "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA"))
c.report()
