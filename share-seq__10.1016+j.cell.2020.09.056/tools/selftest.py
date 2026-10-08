#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import share_seq
from multimodal_page import render
render(share_seq)
print("PASS  both three-barcode modality libraries build")
