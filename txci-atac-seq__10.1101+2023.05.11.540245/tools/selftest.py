#!/usr/bin/env python3
"""Source-boundary checks for the txci-ATAC model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import txci_atac as T
from checks import Check

c = Check()
c.section("txci-ATAC source boundaries")
c("Tn5ME-A source sequence", "".join(s.top for s in T.tn5_a()),
  "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("Tn5ME-B source sequence", "".join(s.top for s in T.tn5_b()),
  "CGTGTGCTCTTCCGATCTGAACCGCGAGATGTGTATAAGAGACAG")
c("Short SBS source sequence", T.SHORT_SBS, "CGTGTGCTCTTCCGATCT")
c("P7.S701 source sequence", "".join(s.top for s in T.p7_primer()),
  "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
c.report()
