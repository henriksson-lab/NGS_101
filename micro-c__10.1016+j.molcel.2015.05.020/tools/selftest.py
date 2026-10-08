#!/usr/bin/env python3
"""Source-boundary checks for the original Micro-C workflow."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import micro_c as M
from checks import Check, run_common

check=Check(); check.section("Hsieh et al. Micro-C protocol choices")
check("both labelled fill-in bases",M.REPAIR.substituted_bases,frozenset(("A","C")))
check("250–350-bp gel selection",M.SIZE_SELECTION_BP,(250,350))
check("12–15 library PCR cycles",M.PCR_CYCLES,(12,15))
check("heterogeneous fill strand is derived",M.REPAIR.fill_5p,"GCTA")
run_common(check); check.report()
