#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telopcr as M
from checks import Check,run_common
c=Check(); c.section("TeloPCR-seq source transcription"); c("S. pombe repeat",M.END.repeat,"GGTTACA"); c("published terminal anchor",M.ANCHOR,"AGGGTTAATAAGCGGCCGCGTCG"); c("published anchor PCR primer",M.ANCHOR_PCR_PRIMER,"CGACGCGGCCGCTTATTAACCCT"); run_common(c); c.report()
