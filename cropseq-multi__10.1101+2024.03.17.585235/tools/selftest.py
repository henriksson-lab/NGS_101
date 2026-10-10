#!/usr/bin/env python3
"""Guide identity and workflow continuity are enforced by the shared constructor."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
chem.workflow()
