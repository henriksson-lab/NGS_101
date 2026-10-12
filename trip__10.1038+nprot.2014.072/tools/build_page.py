#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import trip
from multimodal_page import render

OUT = HERE.parent / "trip.html"

if __name__ == "__main__":
    OUT.write_text(render(trip), encoding="utf-8")
    print(f"wrote {OUT}")
