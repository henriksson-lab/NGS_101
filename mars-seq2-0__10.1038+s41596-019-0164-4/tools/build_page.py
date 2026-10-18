#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem as model
from batch_page import render
OUT = HERE.parent / "mars-seq2-0.html"
if __name__ == "__main__":
    OUT.write_text(render(model), encoding="utf-8")
    print(f"wrote {OUT}")
