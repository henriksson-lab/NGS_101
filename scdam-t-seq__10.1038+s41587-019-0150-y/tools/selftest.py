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
c("printed RPi CGTGAT example", P.RNA_PCR_INDEX,
  "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA")
c("printed RP1", P.RNA_PCR1,
  "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA")
P.protocol_graph().require_publishable()
c.report()
