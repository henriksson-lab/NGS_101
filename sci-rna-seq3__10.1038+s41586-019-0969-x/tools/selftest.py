#!/usr/bin/env python3
"""Source-transcription checks for sci-RNA-seq3 oligos."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import sci_rna_seq3 as S
from checks import Check
c=Check(); c.section("2019 Supplementary Table S11 transcription")
c("sc_ligation_RT_1",S.RT1,"CAGAGCNNNNNNNNTCCTACCAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
c("sc_ligation_1",S.HAIRPIN1,"GCTCTGACAATCAAGTUACGACGCTCTTCCGATCTACTTGATTGT")
c("P5-1",S.P5_1,"AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("P7-1",S.P7_1,"CAAGCAGAAGACGGCATACGAGATCCGAATCCGAGTCTCGTGGGCTCGG")
c.report()
