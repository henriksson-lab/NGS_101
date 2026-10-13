#!/usr/bin/env python3
# Molecular topology, source coverage and sequencing readiness are construction invariants.
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import pna

pna.protocol_graph()
print("PNA model constructs and is publication-complete")
