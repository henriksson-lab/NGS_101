#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "scrb-seq__10.1101+003236" / "tools"
sys.path.insert(0, str(SHARED))
from page_parts import render_page
import scrb

OUT = HERE.parent / "mcscrb-seq.html"
if __name__ == "__main__":
    render_page(scrb.MCSCRB, OUT)

