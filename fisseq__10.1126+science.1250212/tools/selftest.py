#!/usr/bin/env python3
"""Source-boundary checks for the Lee et al. FISSEQ model."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import fisseq as F
from checks import Check, run_common
check = Check()
check.section("Lee et al. protocol choices")
check("RT primer as published", F.RT_PRIMER_WRITTEN,
      "/5phos/TCTCGGGAACGCTGAAGANNNNNN")
check("RCA primer as published", F.RCA_PRIMER_WRITTEN, "TCTTCAGCGTTCCCGA*G*A")
check("CircLigase II incubation", (F.CIRCLIGASE_C, F.CIRCLIGASE_MIN), (60, 60))
run_common(check)
check.report()
