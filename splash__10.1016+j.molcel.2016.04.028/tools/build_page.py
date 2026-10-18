#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path.insert(0,str(H.parents[1]/"lib"))
from protocol_build import build
import chem
OUT=H.parent/"splash.html"
if __name__=="__main__": build(chem,OUT)
