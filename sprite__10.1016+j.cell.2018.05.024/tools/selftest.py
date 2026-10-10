#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sprite_protocol as M
from checks import Check
c=Check(); c("published tag order",[x.name for x in M.TAGS],["DPM tag","Odd tag 1","Even tag","Odd tag 2","Terminal tag"]); c.report()
