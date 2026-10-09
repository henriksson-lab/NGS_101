#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import pairedtag as M
from checks import Check
c=Check(); c.section("PairedTag Supplementary Table 1 transcription")
c("pMENTs", M.PMENTS, "CTGTCTCTTATACACATCT")
c("AdaptorA", M.ADAPTOR_A, "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("Linker-R02", M.LINKER_R02, "CGAATGCTCTGGCCTCTCAAGCACGTGGAT")
c("Linker-R03", M.LINKER_R03, "GGTCTGAGTTCGCACCGAAACATCGGCCAC")
c("published Read 2 barcode windows", M.PUBLISHED_READ2_BARCODE_WINDOWS,
  ((10, 13), (47, 50), (84, 87)))
c.report()
