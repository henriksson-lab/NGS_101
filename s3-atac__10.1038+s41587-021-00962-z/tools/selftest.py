#!/usr/bin/env python3
"""Source-transcription checks for the s3-ATAC protocol-specific oligos."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import s3atac as S
from checks import Check

check = Check()
check.section("Supplementary Tables 2--5 transcription")
check("SBS12_18_UME_sci_1", "".join(s.top for s in S.tn5_transfer()),
      "CGTGTGCTCTTCCGATCTGAACCGCGUAGATGTGTATAAGAGACAG")
check("PCR_i7_P7.S701", "".join(s.top for s in S.i7_pcr_primer()),
      "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
check("PCR_A_i5_A", "".join(s.top for s in S.i5_pcr_primer()),
      "AATGATACGGCGACCACCGAGATCTACACCCTTAAGATCGTCGGCAGCGTC")
check.report()
