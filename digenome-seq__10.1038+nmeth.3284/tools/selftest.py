#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem

c = Check()
# No method-specific WGS oligos are assigned in this capped source scope.
# Construct the native source-limited graph; its delegated library step stays explicit.
chem.protocol_graph()
c.report()
