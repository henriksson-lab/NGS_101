#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import sciatac

OUT = Path(__file__).resolve().parents[1] / "sci-atac-seq.html"
if __name__ == "__main__":
    render_page(sciatac.SCI15, OUT)

