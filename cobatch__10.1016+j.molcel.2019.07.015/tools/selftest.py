#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import cobatch as M
from checks import Check
c=Check()
c("T5 invariant arms", M.T5_OUTER + M.T5_INNER,
  "TCCACGCGCGATCGAGGACGGC")
c("T7 invariant arms", M.T7_OUTER + M.T7_INNER,
  "CTGTCCCTGTCCCACCGTCTCCGCCTC")
c("representative T5 barcodes", M.T5_BARCODES,
  ("TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA"))
c("representative T7 barcodes", M.T7_BARCODES,
  ("CGAGTAAT", "TCTCCGGA", "AATGAGCG", "GGAATCTC"))
c("common ME bottom", M.TN5_ME_REV, "CTGTCTCTTATACACATCT")
c("four structured cell-index parts",
  sum(s.feature is not None and s.feature.role=="cell_barcode"
      for s in M.FINAL_LIBRARY), 4)
c.report()
