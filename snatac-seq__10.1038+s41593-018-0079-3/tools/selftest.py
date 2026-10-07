#!/usr/bin/env python3
"""Source-transcription checks for snATAC-seq Supplementary Table 5."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import snatac as S
from checks import Check

check = Check()
check.section("canonical assembly against representative Supplementary Table 5 rows")
sequence = lambda segments: "".join(s.top for s in segments)
check("assembled p5_1 reproduces the source row",
      sequence(S.p5_transposon(S.P5_1_BARCODE)), S.P5_1)
check("assembled p7_1 reproduces the source row",
      sequence(S.p7_transposon(S.P7_1_BARCODE)), S.P7_1)
check("assembled N701 reproduces the source row",
      sequence(S.p7_index_primer(S.N701_ORDERED_INDEX)), S.N701)
check("assembled S502 reproduces the source row",
      sequence(S.p5_index_primer(S.S502_INDEX)), S.S502)
check.report()
