#!/usr/bin/env python3
"""Run the shared checks for the related SMART-seq constructs."""
from __future__ import annotations

import runpy
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools" / "selftest.py"
runpy.run_path(str(SHARED), run_name="__main__")
