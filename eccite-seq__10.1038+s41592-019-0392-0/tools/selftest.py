#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import eccite_seq
from multimodal_page import render
render(eccite_seq)
print("PASS  transcript, tag and direct-guide libraries build")
