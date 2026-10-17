#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem  # construction validates the final library and primer sites
from checks import Check
c=Check()
# Supplementary Table 3 transcriptions (common tails and one complete locus set).
c("forward primer adaptor tail", chem.FORWARD_TAIL,
  "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("reverse primer adaptor tail", chem.REVERSE_TAIL,
  "TCGGCATTCCTGCTGAACCGCTCTTCCGATCT")
c("MUL1_rs6700034 RT primer", chem.RT_PRIMER, "CTTGAGCTCTTGGGAGACC")
c("MUL1_rs6700034 forward primer", chem.FORWARD_LOCUS, "GAGGATGAACTGGCACAGC")
c("MUL1_rs6700034 reverse primer", chem.REVERSE_LOCUS, "CTGGGCGACTGATTTCGAG")
# PCR amplicons carry no ligation-only dA junction.
c("amplicon has no invented dA junction",
  "dA junction" in [s.name for s in chem.FINAL_LIBRARY], False)
c.report()
