#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check
from chemdraw import revcomp
c=Check(); c("P8 forward target",chem.TARGET_F,"CCAGCAACTGTTTGTGGACCTA")
# Methods Step 4 prints the 102-nt reference used for read matching, barcodes as N.
c("102-nt reference (Methods Step 4)","N"*5+chem.TARGET_AMPLICON+"N"*5,
  "NNNNNACACCTGTGCCTGTTAAACCATTGAAGTTGAAATTGACACATTTGTTTTTAACCAAATTAGTAGACTTTTTAGGTCCACAAACAGTTGCTGGNNNNN")
c("amplicon ends at the forward-primer site",chem.TARGET_AMPLICON.endswith(revcomp(chem.TARGET_F)),True)
c("amplicon starts at the reverse-primer site",chem.TARGET_AMPLICON.startswith(revcomp(chem.TARGET_R)),True)
c.report()
