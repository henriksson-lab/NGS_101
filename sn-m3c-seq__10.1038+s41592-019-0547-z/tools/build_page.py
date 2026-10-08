#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import snm3c as M
from single_cell_hic_page import render
OUT=HERE.parent/"sn-m3c-seq.html"
if __name__ == "__main__": OUT.write_text(render(M), encoding="utf-8"); print(f"wrote {OUT}")
