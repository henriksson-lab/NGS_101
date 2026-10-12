#!/usr/bin/env python3
"""Source-transcription checks for Drop-scChIP-seq."""
from pathlib import Path
import sys

H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
import seqprimers as sp
from checks import Check

c = Check()
c("three A01 top oligos transcribed from Supplementary Table 4",
  tuple(x.top_5to3 for x in P.ROUNDS),
  ("TTCGGCTTGGTTGCGCTTGGTTGC", "TGACACTTGGCGCTACTTGGCGCT",
   "TCCCGCCATTGGACGCCATTGGAC"))
c("common-flank oligos transcribed from Supplementary Table 5",
  (P.SEQ1_TOP, P.SEQ1_BOTTOM, P.SEQ2_TOP, P.SEQ2_BOTTOM),
  ("TTAAGAATTTAATACGACTCACTATAGGGAGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
   "CGAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCTCCCTATAGTGAGTCGTATTAAATTCTTAA",
   "CAACGTGATTGCTTGTGACTTAA", "TTAAGTCACAAGCAATCAC"))
c("final library has three cell identities and no UMI",
  [s.feature.role for s in P.FINAL_LIBRARY if s.feature],
  ["sample_index", "cell_barcode", "cell_barcode", "cell_barcode"])
c("standard run primers bind final library",
  sp.verify(P.FINAL_LIBRARY, P.SEQ_PRIMERS, P.ACQUIRED_ROLES), [])
c.report()
