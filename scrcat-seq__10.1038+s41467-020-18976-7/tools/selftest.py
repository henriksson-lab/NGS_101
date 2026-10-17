#!/usr/bin/env python3
"""Graph construction enforces primer landing, PCR geometry and identifier features."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
chem.protocol_graph()
