#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import nice_library

c = Check()
c.section("NicE-seq paper/kit boundary")
library, primers = nice_library()
c("the unspecified historical kit is visibly inferred",
  all(s.inferred for s in library if s.name != "sonicated genomic DNA carrying an internal biotin-dA or biotin-dC label"), True)
c("no unreported i5 read is invented", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)", "Read 2"))
c.report()
