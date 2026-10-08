#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import directional_rna as M
from checks import Check,run_common
c=Check(); c.section("E7760 defining chemistry"); c("dUTP shown only in second-strand panel",sum("U" in x[0] for r in M.SECOND for x in r.chunks),1); c("all four sequencing reads",len(M.SEQ_PRIMERS),4); run_common(c); c.report()
