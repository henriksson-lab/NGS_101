#!/usr/bin/env python3
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; SHARED=HERE.parents[1]/"split-seq__10.1126+science.aam8999"/"tools"; sys.path.insert(0,str(SHARED))
from page_parts import render_page
OUT=HERE.parent/"microsplit.html"
if __name__=="__main__": OUT.write_text(render_page("microSPLiT"),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
