#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
from checks import Check
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import sci_space_library

c = Check()
c.section("sci-Space molecular boundary")
_, primers = sci_space_library()
c("hybrid sci-RNA-seq run roles", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"))
c.report()
