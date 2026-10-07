#!/usr/bin/env python3
"""Run the shared CEL-Seq construct import validation."""
from __future__ import annotations

import runpy
from pathlib import Path

HERE = Path(__file__).resolve().parent
runpy.run_path(str(HERE.parents[1] / "cel-seq__10.1016+j.celrep.2012.08.003" /
                   "tools" / "selftest.py"), run_name="__main__")
