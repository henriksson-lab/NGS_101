#!/usr/bin/env python3
"""Source-transcription checks for Tang et al. 2009 oligos."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import tang2009 as T
from checks import Check

check = Check()
check.section("2009 supplementary oligo table transcription")
check("UP1", T.UP1, "ATATGGATCCGGCGCGCCGTCGAC" + "T" * 24)
check("UP2", T.UP2, "ATATCTCGAGGGCGCGCCGGATCC" + "T" * 24)
check("P1 top", T.P1_TOP, "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT")
check("P1 bottom", T.P1_BOTTOM, "ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTT")
check("P2 top", T.P2_TOP, "AGAGAATGAGGAACCCGGGGCAGTT")
check("P2 bottom", T.P2_BOTTOM, "CTGCCCCGGGTTCCTCATTCTCT")
check("library PCR 1", T.LIBRARY_PCR_1, "CCACTACGCCTCCGCTTTCCTCTCTATG")
check("library PCR 2", T.LIBRARY_PCR_2, "CTGCCCCGGGTTCCTCATTCT")
check.report()
