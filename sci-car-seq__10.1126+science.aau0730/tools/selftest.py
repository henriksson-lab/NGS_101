#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sci_car as M
from checks import Check
c=Check(); c("two distinct modality libraries",[x[0] for x in M.FINAL_LIBRARIES],["RNA library","ATAC library"]); c.report()
