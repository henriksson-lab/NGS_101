#!/usr/bin/env python3
"""Construction invariants are enforced by the local reaction graph."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem

chem.protocol_graph()
