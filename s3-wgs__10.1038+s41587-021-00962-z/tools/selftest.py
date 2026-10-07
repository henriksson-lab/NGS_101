#!/usr/bin/env python3
"""Structural checks for the s3-WGS construct model."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import s3wgs as S
from checks import Check

check = Check()

check.section("source transcription")
check("representative U-ME oligo matches Supplementary Table 2",
      "".join(s.top for s in S.tn5_transfer_segments()),
      "CGTGTGCTCTTCCGATCTGAACCGCGUAGATGTGTATAAGAGACAG")
check("representative i7 primer matches Supplementary Table 4", S.i7_pcr_primer(),
      "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
check("representative i5 primer matches Supplementary Table 5", S.i5_pcr_primer(),
      "AATGATACGGCGACCACCGAGATCTACACCCTTAAGATCGTCGGCAGCGTC")

check.report()
