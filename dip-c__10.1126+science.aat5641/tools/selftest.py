#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dipc as M
from checks import Check, run_common
c=Check(); c.section("Dip-C protocol choices"); c("drawn enzyme",M.WORKFLOW.enzyme.name,"MboI"); c("no biotin",M.WORKFLOW.marking,"none"); c("amplification",M.WORKFLOW.amplification,"META"); c("META sequence before genomic read",M.META_READ_PREFIX_NT,39); run_common(c); c.report()
