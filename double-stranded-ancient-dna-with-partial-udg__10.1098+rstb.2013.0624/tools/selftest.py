#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import partial_udg as M
from checks import Check
c=Check(); c.section("source transcription"); c("representative P5/P7 barcode pair",(M.P5_BARCODE,M.P7_BARCODE),("ATCGATT","GACTTAT")); c.report()
