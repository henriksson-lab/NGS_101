#!/usr/bin/env python3
"""Source-transcription checks for LAST-seq."""
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import last_seq as P
from checks import Check
c = Check()
c.section("verbatim oligos and stated identifier lengths")
c("RandomRT primer6", P.RANDOM_RT, "TACACGACGCTCTTCCGATCTNNNNNN")
c("dTrU module tail", P.DTRU_TAIL,
  "CCACCTTTCATTCACCCTTTTTTTUUUUUUUUUUUUUUUUUUU")
c("Table S3 NEBNext1 index primer (3' phosphorothioate dropped)",
  P.NEBNEXT_INDEX_5 + "CGTGAT" + P.NEBNEXT_INDEX_3,
  "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
c("Table S3 NEBNext Universal (3' phosphorothioate dropped)", P.NEBNEXT_UNIVERSAL,
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("Table S3 P5 and P7", (P.LIB_P5, P.LIB_P7),
  ("AATGATACGGCGACCACCGAGATCTACAC", "CAAGCAGAAGACGGCATACGAGAT"))
c("cell barcode and UMI lengths", (P.CELL_BARCODE_LEN, P.UMI_LEN), (6, 8))
c.report()
