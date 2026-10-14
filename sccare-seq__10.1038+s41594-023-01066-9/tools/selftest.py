#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sccare_seq as M
from checks import Check
c=Check()
c("all eight RT barcodes transcribed", M.RT_BARCODES,
  ("CTTAGGAC", "GTTGGAGC", "TCTTACGC", "GATGGTCA",
   "GAGAATCC", "CGCAACGA", "CGGACCTA", "GAATTAAC"))
c("representative indexed primers use exact source members",
  (M.N501, M.N701, M.T501),
  ("AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC",
   "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG",
   "AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"))
c.report()
