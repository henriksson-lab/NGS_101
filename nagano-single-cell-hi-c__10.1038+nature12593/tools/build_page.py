#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import nagano as M
from single_cell_hic_page import render
OUT=HERE.parent/"nagano-single-cell-hi-c.html"
if __name__ == "__main__": OUT.write_text(render(M), encoding="utf-8"); print(f"wrote {OUT}")
