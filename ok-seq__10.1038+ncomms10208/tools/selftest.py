#!/usr/bin/env python3
"""Source-transcription checks for OK-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("A1 Watson transcription",P.A1_WATSON,"ACACTCTTTCCCTACACGACGCTCTTCC")
c("A2 Watson transcription",P.A2_WATSON,"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
c("original forward PCR starts at printed base",P.PCR_FORWARD[:8],"ATGATACG")
c.report()
