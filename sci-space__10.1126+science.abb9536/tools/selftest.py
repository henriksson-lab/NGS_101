#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
from spatial_protocols import sci_space_library

c = Check()
c.section("sci-Space molecular boundary")
_, primers = sci_space_library()
c("hybrid sci-RNA-seq run roles", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"))
c.report()
