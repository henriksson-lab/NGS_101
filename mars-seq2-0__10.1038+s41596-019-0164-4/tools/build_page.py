#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "mars-seq__10.1126+science.1247651" / "tools"
sys.path.insert(0, str(SHARED))
from page_parts import render_page
OUT = HERE.parent / "mars-seq2-0.html"
if __name__ == "__main__":
    OUT.write_text(render_page("MARS-seq2.0"), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
