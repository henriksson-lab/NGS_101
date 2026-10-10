#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import epidamid as P
from checks import Check
c=Check();c.section("EpiDamID oligo transcription");c("DamID2 adapter",(P.ADRT,P.ADRB,P.ADR_PCR),("CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA","TCCTCGGCCGCG","GGTCGCGGCCGAGGATC"));c.report()
