#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tea_seq as M
from checks import Check
c=Check();c("three analyte libraries",len(M.FINAL_LIBRARIES),3);c.report()
