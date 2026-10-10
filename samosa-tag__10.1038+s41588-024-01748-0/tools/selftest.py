#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import samosa_tag as M
from checks import Check
c=Check(); c("hairpin Tn5 endpoint is a closed dumbbell", M.SMRTBELL.name, "hairpin-tagmented SMRTbell"); c.report()
