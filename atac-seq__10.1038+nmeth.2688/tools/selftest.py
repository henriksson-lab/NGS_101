#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import atac_seq
from checks import Check
atac_seq.protocol_graph()
Check().report()
