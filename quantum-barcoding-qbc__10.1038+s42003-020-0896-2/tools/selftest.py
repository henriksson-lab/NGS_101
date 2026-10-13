#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import qbc as M
from checks import Check
c=Check(); c.section("supplementary oligo transcription")
c("QBC anchor",M.ANCHOR,"GCTCCCTGTCTGACG")
c("four subcodes plus UMI",(4*7,10),(28,10))
c("antibody core length",len(M.indexed_core(True).top()),133)
M.protocol_graph(); c.report()
