#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import multiome as M
from multimodal_page import render
OUT=H.parent/"multiome.html"
if __name__ == "__main__": OUT.write_text(render(M),encoding="utf-8"); print(f"wrote {OUT}")
