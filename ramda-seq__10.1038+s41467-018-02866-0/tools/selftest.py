#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check
c=Check(); c.section("RamDA-seq source claims")
c("first-strand NSR pool size",408,408)
c("defining study NextSeq read cycles",76,76)
c.report()
