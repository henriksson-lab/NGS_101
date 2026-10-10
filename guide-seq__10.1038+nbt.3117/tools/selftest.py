#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import guide_seq as M
from checks import Check
c=Check();c("molecular barcode is structured",M.FINAL_LIBRARY.get("8-nt molecular barcode").feature.role,"umi");c.report()
