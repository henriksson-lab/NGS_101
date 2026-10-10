#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import m6a_seq as M
from checks import Check
c=Check();c("IP construct is not assigned a molecular barcode",sum(s.feature is not None and s.feature.role not in ("sample_index",) for s in M.FINAL_LIBRARY),0);c.report()
