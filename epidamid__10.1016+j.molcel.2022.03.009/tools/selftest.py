#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import epidamid as P
from checks import Check
c=Check();c.section("EpiDamID oligo transcription")
c("Vogel 2007 m6A-PCR adapter (EpiDamID Key Resources)",(P.ADRT,P.ADRB,P.ADR_PCR),("CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA","TCCTCGGCCGCG","GGTCGCGGCCGAGGATC"))
c("cited DamID2 adapter example (Markodimitraki 2020 Reagents)",(P.DAMID2_TOP,P.DAMID2_BOTTOM),("GGTGATCCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATCNNNTGCANNNTATGGA","TCCATANNNTGCANNNGATCGTCGGACTGTAGAACTCTGAACCCCTATAGTGAGTCGTATTACCGGGAGCTT"))
import build_page
build_page.render(build_page.P)
c.report()
