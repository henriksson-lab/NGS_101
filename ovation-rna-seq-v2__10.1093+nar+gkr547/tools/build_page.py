#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from protocol_build import build
OUT = HERE.parent / "ovation-rna-seq-v2.html"
if __name__ == "__main__": build(__import__("ovation_rna_seq_v2"), OUT)
