#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tail_seq as M
from checks import Check
c=Check();c("tail-facing read is long",M.READ_LENGTHS["Read 2"],251);c.report()
