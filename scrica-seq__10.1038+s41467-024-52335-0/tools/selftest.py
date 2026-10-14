#!/usr/bin/env python3
"""Source-transcription checks for scRICA-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("isoform and accessibility outputs",len(P.FINAL_LIBRARIES),2)
c("linker is the two exact end homologies", P.LINKER,
  "ACTCTGCGTTGATACCACTGCTTGCAATGAAGTCGCAGGGTTG")
c("source-printed TSS enrichment primer", P.I7_ANTI_TAG1,
  "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCGCGTGCCGTATAGCGTTAACT")
c("source-printed terminal primer pools", (P.INDEX_N7, P.INDEX_N5),
  ("CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
   "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"))
c.report()
