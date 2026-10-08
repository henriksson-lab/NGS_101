#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import slide_seqv2
from multimodal_page import render
render(slide_seqv2)
print("PASS  published bead architecture and primer geometry build")
