#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import flex
from multimodal_page import render
render(flex)
print("PASS  adjacent-probe and final-primer geometry build")
