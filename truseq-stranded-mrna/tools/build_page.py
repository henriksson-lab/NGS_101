#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from batch_page import render

OUT = HERE.parent / "truseq-stranded-mrna.html"
OUT.write_text(render(chem), encoding="utf-8")
print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
