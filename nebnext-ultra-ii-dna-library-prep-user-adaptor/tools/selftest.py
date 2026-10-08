#!/usr/bin/env python3
"""Source-boundary checks for the NEBNext USER-adaptor workflow."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import nebnext_user as N
from checks import Check, run_common
check = Check()
check.section("NEB manuals' protocol choices")
check("end-prep incubations", N.END_PREP_TIMES, ((20, 30), (65, 30)))
check("hairpin ligation", N.LIGATION, (20, 15))
check("USER incubation", N.USER, (37, 15))
check("minimum PCR cycles", N.MINIMUM_PCR_CYCLES, 3)
run_common(check)
check.report()
