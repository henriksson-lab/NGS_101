#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import indrop
OUT = Path(__file__).resolve().parents[1] / "indrop-v1.html"
if __name__ == "__main__": render_page(indrop.V1, OUT)

