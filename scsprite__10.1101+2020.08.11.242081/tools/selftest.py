#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scsprite as M
from checks import Check
c=Check(); c("cell barcode has three rounds",len(M.IDENTIFIERS[0].segments),3); c("cluster barcode has six rounds",len(M.IDENTIFIERS[1].segments),6); c.report()
