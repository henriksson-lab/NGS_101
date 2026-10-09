#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import snmc_seq2 as M
from checks import Check
c=Check(); c.section("source transcription"); c("representative P5L_AD002 inline barcode",M.INLINE,"CGATGT"); c.report()
