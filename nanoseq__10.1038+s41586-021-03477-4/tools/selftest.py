#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check
import nanoseq_chem as M
c=Check()
c.section("NanoSeq author-manuscript qPCR oligos")
c("NanoqPCR1", "".join(x.top for x in M.NANOQPCR1), "ACACTCTTTCCCTACACGAC")
c("NanoqPCR2", "".join(x.top for x in M.NANOQPCR2), "GTGACTGGAGTTCAGACGTG")
c.report()
