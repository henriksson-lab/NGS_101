#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from polony_spatial import UMI_PRIMER, TSO
c = Check(); c.section("Pixel-seq source transcription")
c("seven-random-base UMI primer", UMI_PRIMER,
  "CTACACGACGCTCTTCCGATCTNNNNNNNCCGCTAATACGACTCACTATAGGG")
c("template-switching oligo with modifications", TSO,
  "AAGCTGGTATCAACGCAGAGTGAATrG+GrG")
c.report()
