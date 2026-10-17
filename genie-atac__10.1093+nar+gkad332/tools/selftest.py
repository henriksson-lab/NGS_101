#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check
c=Check()
c("MEDS-A transcription", chem.MEDS_A, "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("MEDS-REV transcription", chem.MEDS_REV, "CTGTCTCTTATACACATCT")
c("Tn5A PCR primer", chem.TN5A_PRIMER, "TCGTCGGCAGCGTCAG")
c("rs12269414 biotin linear primer", chem.BIOTIN_LINEAR, "CCAGGAGTTGGAACCTGC")
c("rs12269414 nested forward", chem.NESTED_FORWARD,
  "TCGGCATTCCTGCTGAACCGCTCTTCCGATCTCCTCTTGCCTTGGCTGTG")
c("rs12269414 gDNA-branch reverse oligo tail", chem.GDNA_REVERSE_TAIL,
  "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("PCR amplicon has no invented dA junction",
  "dA junction" in [s.name for s in chem.FINAL_LIBRARY], False)
c.report()
