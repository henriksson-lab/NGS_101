#!/usr/bin/env python3
import runpy
from pathlib import Path
HERE=Path(__file__).resolve().parent
runpy.run_path(str(HERE.parents[1]/"dscatac-seq__10.1038+s41587-019-0147-6"/"tools"/"selftest.py"),run_name="__main__")
