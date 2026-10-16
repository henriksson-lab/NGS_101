#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "lib")]
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import chem
from checks import Check

c = Check()
c.section("Fast-ATAC source transcription")
c("Corces 2016 Online Methods: 76×8×8×76", chem.READ_LENGTHS,
  {"Read 1": 76, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 76})
c("Both cited dual-index primer families present", tuple(p.role for p in chem.SEQ_PRIMERS),
  ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"))
c.report()
