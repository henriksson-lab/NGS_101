#!/usr/bin/env python3
"""Source-transcription checks for Seq-Well oligos."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import seq_well as S
from checks import Check
c=Check(); c.section("2017 Supplementary Table 1 transcription")
c("barcoded bead template",S.BEAD_TEMPLATE,"TTTTTTTAAGCAGTGGTATCAACGCAGAGTACJJJJJJJJJJJJNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
c("template-switching oligo bases",S.TSO,"AAGCAGTGGTATCAACGCAGAGTGAATGGG")
c("SMART PCR primer",S.SMART_PCR,"AAGCAGTGGTATCAACGCAGAGT")
c("P5-SMART hybrid bases",S.P5_HYBRID,"AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")
c("custom Read 1 primer",S.CUSTOM_R1,"GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")
c.report()
