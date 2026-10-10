#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from protocol_build import build
OUT = HERE.parent / "sccc.html"
build(chem, OUT)
