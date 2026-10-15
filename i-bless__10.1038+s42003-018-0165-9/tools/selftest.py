#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import BLESS_PROXIMAL_READ_TAG, BLESS_DISTAL_READ_TAG
c = Check()
c.section("i-BLESS source transcription")
c("proximal/distal analysis tags", (BLESS_PROXIMAL_READ_TAG, BLESS_DISTAL_READ_TAG),
  ("TCGAGGTAGTA", "TCGAGACGACG"))
c.report()
