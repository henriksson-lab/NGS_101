#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import emseq as M
from checks import Check,run_common
c=Check(); c.section("E8015 defining chemistry"); c("unmodified C",M.OUTCOME["C"],"U"); c("5mC protected",M.OUTCOME["5mC"],"C"); c("5hmC protected",M.OUTCOME["5hmC"],"C"); c.section("E8015 p. 3 adaptor trimming sequences vs drawn library")
from chemdraw import revcomp
top="".join(x.top for x in M.final_library())
c("Adaptor Read 1 follows the insert on the top strand", M.TRIM_READ1 in top, True)
c("Adaptor Read 2 follows the insert on the bottom strand", revcomp(M.TRIM_READ2) in top, True)
run_common(c); c.report()
