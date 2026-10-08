#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
from checks import Check
c = Check()
c("v2 places cPPT upstream of the guide cassette",
  [x[0] for x in M.vector_map("v2")].index("cPPT/CTS") <
  [x[0] for x in M.vector_map("v2")].index("hU6-guide-scaffold"))
c.report()
