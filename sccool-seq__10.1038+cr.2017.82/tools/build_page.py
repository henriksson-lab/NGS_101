#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scnome-seq__10.7554+eLife.23203" / "tools"))
from page_parts import render_page
import nomecool
OUT = HERE.parent / "sccool-seq.html"
if __name__ == "__main__": render_page(nomecool.SCCOOL, OUT)

