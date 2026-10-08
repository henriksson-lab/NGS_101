#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import duplex_seq
from checks import Check,run_common
c=Check(); run_common(c); c.report()
