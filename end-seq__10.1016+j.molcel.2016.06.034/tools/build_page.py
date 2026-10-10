#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import end_seq as M
from batch_page import render
OUT=H.parent/"end-seq.html"
if __name__ == "__main__": OUT.write_text(render(M),encoding="utf-8"); print(f"wrote {OUT}")
