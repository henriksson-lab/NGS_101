#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import pbat
from checks import Check
c=Check(); c.section("Miura et al. source transcription")
c("BioPEA2", pbat.BIOPEA2, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("PE-reverse", pbat.PE_REVERSE, "CAAGCAGAAGACGGCATACGAGAT")
c("Primer-3", pbat.PRIMER3,
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c.report()
