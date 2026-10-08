#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dnbseq as M
from batch_page import render
OUT=HERE.parent / "dnbseq.html"
if __name__ == "__main__": OUT.write_text(render(M),encoding="utf-8"); print(f"wrote {OUT}")
