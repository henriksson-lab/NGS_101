#!/usr/bin/env python3
"""Checks for verbatim MATQ-seq oligos."""
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import matq_seq as P
from checks import Check
c = Check(); c.section("supplementary protocol oligos")
c("GAT27", P.GAT27, "GTGAGTGATGGTTGAGGATGTGTGGAG")
c("GAT27dT", P.GAT27_DT,
  "GTGAGTGATGGTTGAGGATGTGTGGAGNNNNNTTTTTTTTTTTTTTTTTTTT")
c("GAT21-6N3G", P.GAT21_6N3G, "GATGGTTGAGGATGTGTGGAGNNNNNNGGG")
c("3NGAT24", P.THREE_N_GAT24, "NNNAGTGATGGTTGAGGATGTGTGGAG")
c.report()
