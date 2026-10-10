#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sccare_seq as M
from checks import Check
c=Check();c("Oligo-dT-1 printed barcode",M.RT1,"CTTAGGAC");c.report()
