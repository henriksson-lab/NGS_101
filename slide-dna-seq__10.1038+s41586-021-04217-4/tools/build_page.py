#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import slide_dna_seq as P
from batch_page import render
OUT=H.parent/"slide-dna-seq.html"
if __name__=="__main__":OUT.write_text(render(P),encoding="utf-8");print(f"wrote {OUT}")
