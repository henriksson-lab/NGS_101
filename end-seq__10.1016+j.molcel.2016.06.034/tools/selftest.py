#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import end_seq as M
from checks import Check
from seqprimers import verify
c=Check(); c("published biotin-dT coordinates are thymines", [M.A1.sequence[i] for i in M.A1.biotin_dt], ["T","T"]); c("final primer sites", verify(M.FINAL_LIBRARY,M.SEQ_PRIMERS), []); c.report()
