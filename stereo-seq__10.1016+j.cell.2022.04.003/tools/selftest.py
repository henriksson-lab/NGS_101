#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import stereo_seq
from batch_page import render
render(stereo_seq)
print("PASS  capture-circle-DNB topology builds")
