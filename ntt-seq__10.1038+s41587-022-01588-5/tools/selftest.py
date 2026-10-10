#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ntt_seq as P
from checks import Check
c=Check();c.section("NTT-seq source transcription");c("MEDSA_1",P.MEDSA1,"TCGTCGGCAGCGTCGGATTGCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG");c("custom primers",(P.CUSTOM_R1,P.CUSTOM_I5),("GCGATCGAGGACGGCAGATGTGTATAAGAGACAG","CTGTCTCTTATACACATCTGCCGTCCTCGATCGC"));c.report()
