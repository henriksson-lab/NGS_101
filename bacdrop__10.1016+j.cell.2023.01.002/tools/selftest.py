#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import bacdrop as M
import seqprimers as sp
from checks import Check
c=Check();c("plate and droplet cell-barcode parts",sum(s.feature is not None and s.feature.role=="cell_barcode" for s in M.FINAL_LIBRARY),2)
c("paper barcode and UMI read windows",
  [(x.read,x.feature.part or x.feature.role,x.cycle_start,x.cycle_end)
   for x in sp.feature_spans(M.FINAL_LIBRARY,M.SEQ_PRIMERS,M.READ_LENGTHS)
   if x.feature.role in ("cell_barcode","umi")],
  [("Read 1","plate well",1,13),("Read 1","umi",14,21),
   ("Index 2 (i5)","droplet",1,16)])
c.report()
