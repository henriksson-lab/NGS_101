#!/usr/bin/env python3
"""Source-boundary checks for SMRTbell prep kit 3.0."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import smrtbell as S
from checks import Check, run_common
check = Check()
check.section("PacBio kit 3.0 protocol choices")
check("repair/end-prep incubations", S.END_PREP_TIMES, ((37, 30), (65, 5)))
check("adapter ligation", S.LIGATION, (20, 30))
check("nuclease cleanup", S.NUCLEASE, (37, 15))
run_common(check)
check.report()
