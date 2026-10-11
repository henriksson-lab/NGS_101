#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check
from illumina_vendor_workflows import get

c = Check()
m = get("truseq-stranded-mrna")
c("dUTP marking precedes selective enrichment",
  [r.action for r in m.workflow().reactions][-2:],
  ["End-repair and add one 3-prime dA", "Ligate TruSeq adapters and selectively PCR-enrich"])
c.report()
