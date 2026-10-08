#!/usr/bin/env python3
"""Source-transcription checks for oligos printed by Smallwood et al."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import scbs as S
from checks import Check

check = Check()
check.section("Online Methods transcription")
check("oligo1", S.OLIGO1_HANDLE + "N" * S.RANDOM_NT,
      "CTACACGACGCTCTTCCGATCTNNNNNNNNN")
check("oligo2", S.OLIGO2_HANDLE + "N" * S.RANDOM_NT,
      "TGCTGAACCGCTCTTCCGATCTNNNNNNNNN")
check("PE1.0", S.PE1,
      "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
check.section("Quail Supplementary Table 1 transcription")
check("dedicated iPCRTag index primer", S.INDEX_PRIMER,
      "AAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC")
check.report()
