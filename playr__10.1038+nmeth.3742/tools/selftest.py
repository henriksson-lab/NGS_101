#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import playr as M
from checks import Check
c=Check(); c.section("published system-1 oligos")
c("insert 1",M.INSERT,"ACGACTGAGTTTGGTCACGAT")
c("universal backbone",M.BACKBONE,"ATTAGCGTCCAGTGAATGCGAGTCCGTCTAGGAGAGTAGTACAGCAGCCGTCAAGAGTGTC")
M.workflow(); c.report()
