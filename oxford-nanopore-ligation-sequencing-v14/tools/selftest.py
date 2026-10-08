#!/usr/bin/env python3
"""Source-boundary checks for SQK-LSK114."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import nanopore_ligation as N
from checks import Check, run_common
check = Check()
check.section("Oxford Nanopore V14 protocol choices")
check("LA top as disclosed", N.LA_TOP, "CCTGTACTTCGTTCAGTTACGTATTGCT")
check("LA bottom as disclosed", N.LA_BOTTOM, "GCAATACGTAACTGAACGAAGTACAGG")
check("repair/end-prep incubations", N.END_PREP_TIMES, ((20, 5), (65, 5)))
run_common(check)
check.report()
