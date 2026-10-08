#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import emseq as M
from checks import Check,run_common
c=Check(); c.section("E8015 defining chemistry"); c("unmodified C",M.OUTCOME["C"],"U"); c("5mC protected",M.OUTCOME["5mC"],"C"); c("5hmC protected",M.OUTCOME["5hmC"],"C"); run_common(c); c.report()
