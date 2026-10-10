#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import music as M
from checks import Check
c=Check();c("RNA-linker fixed prefix",M.RNA_LINKER[:12],"CGAGGAGCGCTT");c("three cell rounds",sum(1 for s in M.FINAL if s.feature and s.feature.role=="cell_barcode"),3);c.report()
