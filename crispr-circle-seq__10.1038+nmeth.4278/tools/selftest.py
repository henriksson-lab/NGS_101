#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import crispr_circle as M
from checks import Check
c=Check();c("published stem-loop contains USER site",M.STEM_LOOP.count("U"),1);c.report()
