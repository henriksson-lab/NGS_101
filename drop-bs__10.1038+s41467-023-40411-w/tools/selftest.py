#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import drop_bs as P
from checks import Check
c=Check();c.section("Drop-BS source transcription");c("N9 random primer",P.RANDOM_PRIMER,"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNN");c("cell barcode length",P.CELL_LEN,15);c.report()
