#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scnome-seq__10.7554+eLife.23203" / "tools"))
from selftest_common import run
import nomecool
if __name__ == "__main__": run(nomecool.SCCOOL)

