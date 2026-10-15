#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check
import chem

c = Check()
m = chem
actions = [r.action for r in m.workflow().reactions]
c("dUTP marking precedes selective enrichment",
  actions.index("RNase H / DNA polymerase I second-strand synthesis with dUTP")
  < actions.index("Ligate TruSeq adapters and selectively PCR-enrich"),
  True)
c.report()
