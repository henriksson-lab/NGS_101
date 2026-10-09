#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ddrad as D
from checks import Check
c=Check(); c.section("source transcription")
c("published P1.1", D.P1_TOP, "ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATG")
c("published P1.2", D.P1_BOTTOM, "AATTCATGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
c("published P2.1", D.P2_TOP, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
c("published P2.2", D.P2_BOTTOM, "CGAGATCGGAAGAGCGAGAACAA")
c.report()
