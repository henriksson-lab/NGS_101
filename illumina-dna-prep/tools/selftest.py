#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check
import chem

c = Check()
m = chem
c("bead-linked tagmentation retains its wash",
  [r.action for r in m.workflow().reactions][1:3],
  ["On-bead tagmentation", "Post-tagmentation wash"])
c.report()
