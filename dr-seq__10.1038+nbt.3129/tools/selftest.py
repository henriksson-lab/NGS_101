#!/usr/bin/env python3
"""Source-transcription checks for DR-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check

c = Check()
c("Ad-1x fixed prefix", P.AD1_SPACER + P.T7_PROMOTER + P.RA5,
  "CGATTGAGGCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATC")
c("Ad-1x cell barcode", P.CELL_BARCODE_NT, 8)
c("Ad-1x dT", P.POLYT_NT, 24)
c("Ad-2 handle", P.AD2_HANDLE, "GTGAGTGATGGTTGAGGTAGTGTGGAG")
c("Ad-2 random foot", P.AD2_RANDOM_NT, 8)
c("P1", P.P1, "CGATTGAGGCCGGTAATAC")
c("P2", P.P2, "GTGAGTGATGGTTGAGGTAGTGTGGAG")
c("biotin-P3 bases", P.P3, "GTGAGCTGGAGTTGAGGTAGTGTGGAG")
c("genome and transcriptome outputs", len(P.FINAL_LIBRARIES), 2)
c("RNA identity is a cell barcode, not a UMI",
  [s.feature.role for s in P.RNA if s.feature],
  ["sample_index", "cell_barcode", "sample_index"])
c.report()
