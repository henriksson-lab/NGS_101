#!/usr/bin/env python3
"""Source-transcription checks for CELLO-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("phosphorylated oligo(dT) transcription", P.DT_OLIGO,
  "TGGCGTAGCGGGTTCGAGCGCACCGCAGGGTATCCGGCTATTTTTTTTTTTTTTTTTTTTTTTTTVN")  # Table S5 oligo_dT verbatim
c("top splint transcription", P.TOP_SPLINT,
  "TAGCCGGATACCCTGCGGTGCGCTCGAACCCGCTACGCCATGTAAAACGACGGCCAGT")
c("bottom splint idx1 transcription", P.BOTTOM_SPLINT,
  "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
  "AAGAAAGTTGTCGGTGTCTTTGTGAGTGGTATC" + "NRY" * 7 + "N"
  "ACTGGCCGTCGTTTTACATGGCGTAGCGGGTTCGAGCGCACCGCAGGGTATCCGGCTATTTTTTTTTTTTTTT")  # Table S5 idx1 verbatim
c("patterned UMI length",
  next(len(s.top) for s in P.FINAL_LIBRARY if s.feature and s.feature.role=="umi"), 22)
c("cell barcode length", len(P.BARCODE_1), 24)
c.report()
