#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check
import ovation_rna_seq_v2 as model
c=Check(); c("finished graph is publishable", model.protocol_graph().publishable, True); c.report()
