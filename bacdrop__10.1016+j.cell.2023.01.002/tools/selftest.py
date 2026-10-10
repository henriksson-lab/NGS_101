#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import bacdrop as M
from checks import Check
c=Check();c("plate and droplet cell-barcode parts",sum(s.feature is not None and s.feature.role=="cell_barcode" for s in M.FINAL_LIBRARY),2);c.report()
