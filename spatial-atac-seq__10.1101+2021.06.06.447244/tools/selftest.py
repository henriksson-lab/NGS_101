#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check
c=Check(); c("custom Read 1 transcription",chem.CUSTOM_R1,"GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"); c.report()
