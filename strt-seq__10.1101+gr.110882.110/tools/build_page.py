#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import strt

OUT = Path(__file__).resolve().parents[1] / "strt-seq.html"
if __name__ == "__main__":
    render_page(strt.STRT, OUT)

