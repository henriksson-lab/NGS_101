#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sci_plex as P
from checks import Check
c=Check();c.section("sci-Plex identity roles");c("cell identity is structured",any(s.feature and s.feature.role=="cell_barcode" for s in P.FINAL_LIBRARY));c.report()
