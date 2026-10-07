#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "indrop-v1__10.1016+j.cell.2015.04.044" / "tools"))
from selftest_common import run
import indrop
if __name__ == "__main__": run(indrop.V2)

