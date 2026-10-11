#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from protein_dna_mapping import organic_library

c = Check()
c.section("ORGANIC paired-footprint library")
library, primers = organic_library()
c("the source-defined insert is a native MNase footprint",
  next(s.name for s in library if "footprint" in s.name),
  "native MNase footprint retained with the affinity-purified factor")
c("the historical library is paired-end single-index",
  tuple(p.role for p in primers), ("Read 1", "Index 1 (i7)", "Read 2"))
c.report()
