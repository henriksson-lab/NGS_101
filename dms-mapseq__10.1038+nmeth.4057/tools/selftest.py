#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import dms_mapseq as M
from checks import Check
c=Check(); c("published linker-2",M.LINKER,"CACTCGGGCACCAAGGA"); c("cDNA circle closes RT handle to copied RNA",(M.CIRCLE.closure_left,M.CIRCLE.closure_right),("mutation-bearing cDNA","RT-primer handle")); c.report()
