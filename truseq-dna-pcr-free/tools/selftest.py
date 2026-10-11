#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from checks import Check
from illumina_vendor_workflows import get

c = Check()
# Construction validates the canonical library and all declared primer landing sites.
m = get("truseq-dna-pcr-free")
c("PCR-free workflow has no PCR reaction", any("PCR" in r.action for r in m.workflow().reactions), False)
c.report()
