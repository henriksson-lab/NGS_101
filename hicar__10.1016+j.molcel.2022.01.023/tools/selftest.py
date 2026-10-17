#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from chem import HICAR7_I0_WRITTEN

c = Check()
c.section("HiCAR later-2023 source transcription (2022 assignment unresolved)")
c("HiCAR7-i0 ordering string", HICAR7_I0_WRITTEN,
  "/5Phos/taGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGGAGATGTGTATAAGAGACAG")
c.report()
