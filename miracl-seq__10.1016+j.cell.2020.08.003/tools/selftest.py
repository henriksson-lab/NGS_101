#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H.parents[1]/"lib")]
from checks import Check
from chem import MiraclLoading
c=Check(); c.section("MIRACL construction boundary")
try:
    MiraclLoading(1)
    rejected=False
except ValueError:
    rejected=True
c("ordinary loading is rejected by the MIRACL constructor", rejected, True)
c.report()
