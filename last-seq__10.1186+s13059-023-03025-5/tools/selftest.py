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
c("cell barcode and UMI lengths", (P.CELL_BARCODE_LEN, P.UMI_LEN), (6, 8))
c.report()
