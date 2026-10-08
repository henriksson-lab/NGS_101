#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import snm3c as M
from checks import Check, run_common
c=Check(); c.section("sn-m3C-seq protocol choices"); c("amplification",M.WORKFLOW.amplification,"bisulfite-PCR"); c("no biotin",M.WORKFLOW.marking,"none"); c("random-primer trim",M.RANDOM_PRIMER_TRIM_NT,25); c("Adaptase tail trim",M.ADAPTASE_TAIL_TRIM_NT,3); run_common(c); c.report()
