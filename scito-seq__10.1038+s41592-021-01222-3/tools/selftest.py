#!/usr/bin/env python3
"""Construct the protocol-owned declaration and its graph invariants."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
chem.workflow()

chem.protocol_graph()
