#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path.insert(0,str(H.parents[1]/"lib"))
from multimodal_page import render
import chem
OUT=H.parent/"redeem.html"
if __name__=="__main__": OUT.write_text(render(chem),encoding="utf-8");print(f"wrote {OUT}")
