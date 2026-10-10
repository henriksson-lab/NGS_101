#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import spatial_cut_tag as P
from checks import Check
c=Check();c.section("Spatial-CUT&Tag identifiers");c("two spatial barcode parts",[(s.feature.group,s.feature.part) for s in P.FINAL_LIBRARY if s.feature and s.feature.role=="spatial_barcode"],[("spatial_pixel","row"),("spatial_pixel","column")]);c.report()
