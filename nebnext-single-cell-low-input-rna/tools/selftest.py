#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import low_input_rna as M
from checks import Check,run_common
c=Check(); c.section("E6420 disclosure boundary"); c("RT and TSO handles remain inferred",all(s.inferred for s in M.WTA if "handle" in s.name),True); c("all four sequencing reads",len(M.SEQ_PRIMERS),4); run_common(c); c.report()
