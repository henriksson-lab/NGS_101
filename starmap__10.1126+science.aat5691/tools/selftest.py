#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import starmap as M
from checks import Check
c=Check(); c.section("defining-paper transcription")
c("target-code width", len(M.GENE_CODE), 5)
c("reading-probe width", len(M.DOCK), 11)
M.protocol_graph(); c.report()
