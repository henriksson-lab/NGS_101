#!/usr/bin/env python3
# Construction checks run in chem.protocol_graph(); this suite only confirms it is publishable.
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem
c = Check(); c("finished graph is publishable", chem.protocol_graph().publishable, True); c.report()
