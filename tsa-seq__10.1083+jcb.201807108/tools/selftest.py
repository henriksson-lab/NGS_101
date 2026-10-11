#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
import seqprimers as sp
from checks import Check
from spatial_proximity_labeling import tsa_library

c = Check()
c.section("TSA-seq defining run boundary")
_, primers = tsa_library()
c("single-end, single-index run", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)"))
c.report()
