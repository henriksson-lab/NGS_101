#!/usr/bin/env python3
"""Construction checks run by building the graph; identifier provenance is enforced
by the constructors and primer landing by the sequencing section."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
chem.workflow()
chem.protocol_graph()
