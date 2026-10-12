#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import probac_seq as P
from checks import Check
c=Check()
c.section("ProBac-seq source transcription")
c("extension handle", P.EXTENSION_HANDLE, "AAGCTTTGCTCGAAGAAACTATG")
c("PCR handle", P.PCR_HANDLE, "GTATGGACCGGGTCATC")
c("RCA-1 / RCA-3 primer", P.RCA1_PRIMER, "TTCGAGCAAAGCTTGATGACCCGGTCC")
c("RCA-2 primer", P.RCA2_PRIMER, "GGACCGGGTCATCAAGCTTTGCTCGAA")
c("blocked extension template", P.EXTENSION_TEMPLATE,
  "T"*30+"V"+"N"*11+"AAGCTTGCTCGAAGAAACTATG")
c("10x forward primer", P.PARTIAL_R1, "CTACACGACGCTCTTCCGATCT")
c("InDrop reverse primer", P.INDROP_REVERSE, "GATGACCCGGTCCATACA")
c.report()
