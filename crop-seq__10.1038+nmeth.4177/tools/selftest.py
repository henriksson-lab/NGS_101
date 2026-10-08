#!/usr/bin/env python3
"""Primary-source transcription checks for CROP-seq."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import crop_seq as C
from checks import Check
check = Check()
check.section("CROP-seq primary-source transcription")
check("pooled guide oligo", "".join(s.top for s in C.pooled_guide_oligo()),
      "TGGAAAGGACGAAACACCG" + "N" * 20 + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC")
check("custom Read 1 primer", C.CUSTOM_R1,
      "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")
check.report()
