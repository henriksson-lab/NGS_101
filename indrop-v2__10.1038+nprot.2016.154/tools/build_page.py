#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "indrop-v1__10.1016+j.cell.2015.04.044" / "tools"))
from page_parts import render_page
import indrop
OUT = HERE.parent / "indrop-v2.html"
if __name__ == "__main__": render_page(indrop.V2, OUT)

