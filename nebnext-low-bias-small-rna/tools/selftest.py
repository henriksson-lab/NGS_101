#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import low_bias_small_rna as M
from checks import Check,run_common
c=Check(); c.section("E3420 defining chemistry"); c("compatible 5-prime end","5′-phosphorylated" in next(iter(M.INPUT)).name,True); c("no DNA dA junction","dA junction" in tuple(s.name for s in M.final_library()),False); run_common(c); c.report()
