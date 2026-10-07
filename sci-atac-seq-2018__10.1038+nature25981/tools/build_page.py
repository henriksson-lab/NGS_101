#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "sci-atac-seq__10.1126+science.aab1601" / "tools"
sys.path.insert(0, str(SHARED))
from page_parts import render_page
import sciatac

OUT = HERE.parent / "sci-atac-seq-2018.html"
if __name__ == "__main__":
    render_page(sciatac.SCI18, OUT)

