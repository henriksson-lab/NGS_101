#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import nano_cut_tag as P
from batch_page import render
OUT=H.parent/"nano-cut-tag.html"
if __name__=="__main__":OUT.write_text(render(P),encoding="utf-8");print(f"wrote {OUT}")
