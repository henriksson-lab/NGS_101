#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import emseq as M
from neb_commercial_page import render
OUT=HERE.parent/"nebnext-em-seq-v2.html"
if __name__=="__main__": OUT.write_text(render(M),encoding="utf-8"); print(f"wrote {OUT}")
