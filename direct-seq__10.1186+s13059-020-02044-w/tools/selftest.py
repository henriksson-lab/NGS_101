#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import direct_seq as P
from checks import Check
c=Check();c.section("Direct-seq source transcription")
c("8A8G tract",P.CAPTURE,"AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA")
c("human Gln tRNA", P.HUMAN_TRNA_GLN,
  "GGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCT")
c("nested-PCR reverse primer", P.TRNA_READ2,
  "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTCTGGACTCTGAATCCAGCGAT")
c("guide library preserves one guide identifier",
  [s.feature.role for s in P.GUIDE_LIB if s.feature],
  ["cell_barcode", "umi", "guide_barcode", "sample_index"])
c.report()
