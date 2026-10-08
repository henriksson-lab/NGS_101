#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import scihic as M
from checks import Check, run_common
c=Check(); c.section("sci-Hi-C protocol choices"); c("bridge",M.WORKFLOW.marking,"bridge-adaptor"); c("barcode rounds",M.WORKFLOW.barcoding_rounds,2); c("internal barcode count",sum("round-1 barcode" in s.name for s in M.contact_product()),2); run_common(c); c.report()
