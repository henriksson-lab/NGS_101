#!/usr/bin/env python3
"""Source-transcription checks for Multi-CUT&Tag Table S1."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
import seqprimers

c = Check()
c("A-family antibody barcodes", P.BARCODES_A,
  ("TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA"))
c("B-family antibody barcodes", P.BARCODES_B,
  ("CGAGTAAT", "TCTCCGGA", "AATGAGCG", "GGAATCTC"))
c("Tn5MErev", P.TN5_ME_REV, "CTGTCTCTTATACACATCT")
c("custom Read1", P.READ1, "TCGTCGGCAGCGTCTCCACGC")
c("custom Read2", P.READ2, "GTCTCGTGGGCTCGGCTGTCCCTGTCC")
c("custom Index1", P.INDEX1, "GGACAGGGACAGCCGAGCCCACGAGAC")
c("custom Index2", P.INDEX2, "GCGTGGAGACGCTGCCGACGA")
c("all custom primers land", seqprimers.verify(P.FINAL_LIBRARY, P.SEQ_PRIMERS), [])
c.report()
