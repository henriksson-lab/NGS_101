#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import fs_pcrfree as M
from checks import Check,run_common
c=Check(); c.section("E7430 defining chemistry"); c("UMI length",M.UMI_NT,12); c("i7 plus UMI index read",M.I7_READ_NT,20); c("350-bp fragmentation",M.FS_PROGRAMS[350],(10,37,30,65)); c("450-bp fragmentation",M.FS_PROGRAMS[450],(8,37,30,65)); run_common(c); c.report()
