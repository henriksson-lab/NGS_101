#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chip_seq
from checks import Check
c = Check()
c("Genomic PCR primer 2.1 transcription", chip_seq.GENOMIC_PCR2,
  "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT")
c.report()
