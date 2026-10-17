#!/usr/bin/env python3
"""Construction checks: sensor pairing, primer landing and graph continuity are enforced on build."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from batch_page import render
render(chem)
print("PASS  smol-seq sensor pairing, round-2 primer landing and sequencing-primer geometry")
