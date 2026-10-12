#!/usr/bin/env python3
"""Source-transcription checks for scDam&T-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("printed DamID top example", P.DAM_TOP,
  "GGTGATCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATCNNNTGCANNNTATGGA")
c("printed DamID bottom example", P.DAM_BOTTOM,
  "TCCATANNNTGCANNNGATCGTCGGACTGTAGAACTCTGAACCCCTATAGTGAGTCGTATTACCGGGAGCTT")
c("printed CEL-seq2 example", P.CEL_PRIMER,
  "GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATCNNNGATGNNNTCAT" + "T" * 25 + "V")
c("RNA and DamID outputs",len(P.FINAL_LIBRARIES),2)
c.report()
