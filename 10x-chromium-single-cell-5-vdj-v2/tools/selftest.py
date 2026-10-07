#!/usr/bin/env python3
"""Import-time construct validation for 10x 5' V(D)J v2."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"10x-chromium-single-cell-5-vdj"/"tools"),str(HERE.parents[1]/"lib")]
import vdj_v2  # noqa: F401
from checks import Check
Check().report()
