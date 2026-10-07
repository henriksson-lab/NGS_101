#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
HERE=Path(__file__).resolve().parent; OUT=HERE.parent/"split-seq.html"
if __name__=="__main__": OUT.write_text(render_page("SPLiT-seq"),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
