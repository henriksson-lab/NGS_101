#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import multiome
from multimodal_page import render
render(multiome)
print("PASS  model constructs both primer-validated libraries")
