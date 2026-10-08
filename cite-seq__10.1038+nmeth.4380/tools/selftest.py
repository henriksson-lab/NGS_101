#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import cite_seq
from multimodal_page import render
render(cite_seq)
print("PASS  transcript and ADT primer geometry build")
