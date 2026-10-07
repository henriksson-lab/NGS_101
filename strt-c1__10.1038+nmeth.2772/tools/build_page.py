#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "strt-seq__10.1101+gr.110882.110" / "tools"
sys.path.insert(0, str(SHARED))
from page_parts import render_page
import strt

OUT = HERE.parent / "strt-c1.html"
if __name__ == "__main__":
    render_page(strt.C1, OUT)

