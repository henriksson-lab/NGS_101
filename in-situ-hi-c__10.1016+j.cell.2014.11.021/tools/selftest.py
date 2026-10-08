#!/usr/bin/env python3
"""Source-boundary checks for the Rao et al. in situ Hi-C model."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import insitu_hic as H
from checks import Check, run_common
check = Check()
check.section("Rao et al. protocol choices")
check("restriction enzyme", H.DEMO.enzyme.name, "MboI")
check("biotinylated fill-in nucleotide", H.FILLED.biotin_base, "A")
check("sheared library range", H.SHEAR_RANGE_BP, (300, 500))
run_common(check)
check.report()
