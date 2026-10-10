#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check
c=Check(); c("P8 forward target",chem.TARGET_F,"CCAGCAACTGTTTGTGGACCTA"); c.report()
