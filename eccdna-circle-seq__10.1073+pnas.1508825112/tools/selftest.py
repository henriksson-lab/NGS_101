#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import eccdna_circle as M
from checks import Check
c=Check();c("one RCA unit equals the input circle",len(M.RCA.unit),len(M.CIRCLE));c.report()
