#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import strand_seq
from checks import Check
c = Check()
c("printed Strand-seq index primer",
  strand_seq.INDEX_READ, "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG")
c.report()
