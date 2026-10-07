#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "sci-atac-seq__10.1126+science.aab1601" / "tools"))
from selftest_common import run
import sciatac
if __name__ == "__main__": run(sciatac.SCI3)

