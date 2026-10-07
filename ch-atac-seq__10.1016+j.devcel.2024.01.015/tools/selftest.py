#!/usr/bin/env python3
"""Checks only the external sequence boundaries used by the CH-ATAC model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import ch_atac as C
from checks import Check

c = Check()
c.section("CH-ATAC source boundaries")
c("CH-RNA HY barcode 1", C.HY_BARCODE, "TTCTCGCATG")
c("CH-RNA MGI P7 index 1", C.I7_OLIGO_INDEX, "TAGGTCCGAT")
c("MGI P5", C.MGI_P5, "GAACGACATGGCTACGATCCGACTT")
c("MGI P7", C.MGI_P7, "TGTGAGCCAAGGAGTTGTTGTCTTC")
c.report()
