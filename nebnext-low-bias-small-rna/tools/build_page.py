#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import low_bias_small_rna as M
from neb_commercial_page import render
OUT=HERE.parent/"nebnext-low-bias-small-rna.html"
if __name__=="__main__": OUT.write_text(render(M),encoding="utf-8"); print(f"wrote {OUT}")
