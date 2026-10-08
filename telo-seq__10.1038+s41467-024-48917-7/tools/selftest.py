#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import teloseq as M
from checks import Check,run_common
c=Check(); c.section("Telo-seq source choices"); c("six telorettes",len(M.TELORETTES),6); c("shared core",[x.startswith(M.CORE) for x in M.TELORETTES],[True]*6); c("S1",M.S1,"ACTTCGTTCAGTTACGTATTGCTAGCAAT"); run_common(c); c.report()
