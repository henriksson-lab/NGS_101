#!/usr/bin/env python3
import runpy
from pathlib import Path
HERE = Path(__file__).resolve().parent
runpy.run_path(str(HERE.parents[1] / "mars-seq__10.1126+science.1247651" /
                   "tools" / "selftest.py"), run_name="__main__")
