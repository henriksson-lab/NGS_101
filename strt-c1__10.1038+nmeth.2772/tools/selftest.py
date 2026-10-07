#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "strt-seq__10.1101+gr.110882.110" / "tools"
sys.path.insert(0, str(SHARED))
from selftest_common import run
import strt

if __name__ == "__main__":
    run(strt.C1)

