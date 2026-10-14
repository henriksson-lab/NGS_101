#!/usr/bin/env python3
"""Source-transcription checks for CRISPR-sciATAC.

Every assembled oligo below is compared against one row copied verbatim from
Liscovitch-Brauer et al. 2021, Supplementary Table 2. Everything the model can derive
(complements, decompositions, primer landing sites, read spans) is enforced in
tools/chem.py instead.
"""
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
from checks import Check

c = Check()

# Supplementary Table 2, "Common oligos"
c("T5 tagment sciATAC",
  P.T5_TAGMENT,
  "TCGTCGGCAGCGTCTCCACGCAACCAACGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG")
c("P5 sgRNA PCR2",
  P.P5_SGRNA_PCR2,
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGCCGTCCAGCTGA")
c("P5 sgRNA PCR2 +8 stagger",
  P.P5_SGRNA_PCR2_STAGGERED,
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNCGCCGTCCAGCTGA")

# One row from each barcode plate, with its 8-nt window written back in.
c("Barcoded T7 tagment ATAC, row 1",
  P.T7_TAGMENT.replace(P.BC, "TCTCGCGT"),
  "GTCTCGTGGGCTCGGCTGTCCCTGTCCTCTCGCGTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG")
c("RTi sgRNA, row 1 (same BC1 as MEDS B row 1)",
  P.RT_PRIMER.replace(P.BC, "TCTCGCGT"),
  "CCGGGAGCTGCATGTTCTCGCGTTTTTTCAAGTTGATAA")
c("P5i ATAC PCR1, row 1",
  P.P5I_ATAC.replace(P.BC, "GGTATTGG"),
  "AATGATACGGCGACCACCGAGATCTACACGGTATTGGTCGTCGGCAGCGTC")
c("P5i sgRNA PCR1, row 1 (same BC2 as P5i ATAC row 1)",
  P.P5I_SGRNA.replace(P.BC, "GGTATTGG"),
  "CGCCGTCCAGCTGAGGTATTGGCCGGGAGCTGCATGT")
c("P7i ATAC PCR1, row 1",
  P.P7I_ATAC.replace(P.BC, "AGTAGGTA"),
  "CAAGCAGAAGACGGCATACGAGATAGTAGGTAGTCTCGTGGGCTCGG")
c("P7i sgRNA PCR2, row 1 (same BC3 as P7i ATAC row 1)",
  P.P7I_SGRNA.replace(P.BC, "AGTAGGTA"),
  "CAAGCAGAAGACGGCATACGAGATAGTAGGTATTGTGGAAAGGACGAAACAC")

# Supplementary Fig. 1d,e: the run the paper actually performed.
c("ATAC run lengths", P.ATAC_READ_LENGTHS,
  {"Read 1": 45, "Index 1 (i7)": 43, "Index 2 (i5)": 8, "Read 2": 45})
c("guide run is one unindexed read", P.GUIDE_READ_LENGTHS, {"Read 1": 162})
c.report()
