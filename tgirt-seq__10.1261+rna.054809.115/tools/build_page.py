#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from batch_page import render
import tgirt_seq
OUT=H.parent/"tgirt-seq.html"
if __name__=="__main__": OUT.write_text(render(tgirt_seq),encoding="utf-8"); print(f"wrote {OUT}")
