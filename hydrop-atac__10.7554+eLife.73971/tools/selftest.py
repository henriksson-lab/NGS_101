#!/usr/bin/env python3
"""Source-transcription checks for HyDrop-ATAC protocol-specific oligos."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import hydrop_atac as H
from checks import Check

check = Check()
check.section("Supplementary files 1 and 4 transcription")
check("Acrydite_primer", H.ACRYDITE_PRIMER,
      "TTTTTTTTAATACGACTCACTATAGGGAAGCAGTGGTATCAACGCAGAGTAC")
check("plate 1, first well", H.PLATE1_FIRST, "GCAGTAGCTGTGTAGCAAGTGTACTCTGCG")
check("plate 2, first well", H.PLATE2_FIRST, "AGGGTACTCGTTAGTTGGACGCAGTAGCTG")
check("ATAC plate 3, first well", H.PLATE3_FIRST,
      "CCGAGCCCACGAGACTGACCGTACTAGGGTACTCG")
check("HYi7_1_CGCTCAGTTC", "".join(s.top for s in H.hyi7_primer()),
      "CAAGCAGAAGACGGCATACGAGATCGCTCAGTTCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")
check("HYi5_1_TCGTGGAGCG", "".join(s.top for s in H.hyi5_primer()),
      "AATGATACGGCGACCACCGAGATCTACACTCGTGGAGCGTCGTCGGCAGCGTCAGATGTG")
check.report()
