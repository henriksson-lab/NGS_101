#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check
from illumina_vendor_workflows import get

c = Check()
m = get("nextera-xt-dna-library-prep")
c("solution tagmentation is followed by neutralization",
  [r.action for r in m.workflow().reactions[:2]],
  ["Solution-transposome tagmentation", "Neutralize transposase"])
c.report()
