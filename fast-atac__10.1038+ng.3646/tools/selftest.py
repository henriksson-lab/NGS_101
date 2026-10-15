#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protocol_families.atac import library
from checks import Check

c = Check()
c.section("Fast-ATAC construct boundary")
final, primers = library("Fast-ATAC")
c("historical run roles", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)", "Read 2"))
c.report()
