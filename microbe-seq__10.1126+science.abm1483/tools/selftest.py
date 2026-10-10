#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import microbe_seq as M
from checks import Check
c=Check();c("one structured microbe barcode",sum(s.feature is not None and s.feature.role=="cell_barcode" for s in M.FINAL_LIBRARY),1);c.report()
