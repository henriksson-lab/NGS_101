#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scget_seq as P
from checks import Check
c=Check();c.section("scGET-seq source boundary");c("two identifiers stay distinct",[s.feature.role for s in P.FINAL_LIBRARY if s.feature],["cell_barcode","feature_barcode","sample_index"]);c.report()
