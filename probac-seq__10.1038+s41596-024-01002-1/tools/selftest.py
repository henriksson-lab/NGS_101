#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import probac_seq as P
from checks import Check
c=Check()
c.section("ProBac-seq source transcription")
c("extension handle", P.EXTENSION_HANDLE, "AAGCTTTGCTCGAAGAAACTATG")
c("PCR handle (common 3-prime flank of every Supplementary Table 1 probe)", P.PCR_HANDLE, "TGTATGGACCGGGTCATC")
c("RCA-1 / RCA-3 primer", P.RCA1_PRIMER, "TTCGAGCAAAGCTTGATGACCCGGTCC")
c("RCA-2 primer", P.RCA2_PRIMER, "GGACCGGGTCATCAAGCTTTGCTCGAA")
c("blocked extension template", P.EXTENSION_TEMPLATE,
  "T"*30+"V"+"N"*11+"AAGCTTGCTCGAAGAAACTATG")
c("10x forward primer", P.PARTIAL_R1, "CTACACGACGCTCTTCCGATCT")
c("InDrop reverse primer", P.INDROP_REVERSE, "GATGACCCGGTCCATACA")
c("gel-bead oligo after partial Read 1 (CG000204 Rev D appendix: N16-N12-T30)",
  [(s.name, s.top) for s in P.bead_oligo()[1:]],
  [("cell barcode", "B"*16), ("10x UMI", "U"*12), ("poly(dT)30", "T"*30)])
c.report()
