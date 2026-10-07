#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import nomecool
OUT = Path(__file__).resolve().parents[1] / "scnome-seq.html"
if __name__ == "__main__": render_page(nomecool.SCNOME, OUT)

