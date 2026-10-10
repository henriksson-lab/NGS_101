#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import asap_seq as M
from checks import Check
c=Check(); c("two separately sequenced modalities",len(M.FINAL_LIBRARIES),2); c.report()
