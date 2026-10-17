#!/usr/bin/env python3
"""Source-transcription checks for UDA-seq (Supplementary Table 7, primer summary)."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check

def joined(parts):
    return "-".join("[8-bp sample index]" if s.placeholder else s.top for s in parts)

c=Check()
c.section("UDA-seq Supplementary Table 7 transcription")
c("Truseq-i5 index primer template", joined(P.truseq_i5_index_primer()),
  "AATGATACGGCGACCACCGAGATCT-[8-bp sample index]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("Nextare-i7 index primer template", joined(P.nextera_i7_index_primer()),
  "CAAGCAGAAGACGGCATACGAGAT-[8-bp sample index]-GTCTCGTGGGCTCGG")
c("P5 primer", P.P5_PRIMER, "AATGATACGGCGACCACCGAGA")
c("Partial TSO", P.PARTIAL_TSO, "AAGCAGTGGTATCAACGCAGAG")
c("Tn5-top_ME bases", "".join(s.top for s in P.tn5_top()), "CTGTCTCTTATACACATCT")
c("Tn5-bottom_Read2N", "".join(s.top for s in P.tn5_bottom()), "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")
c.report()
