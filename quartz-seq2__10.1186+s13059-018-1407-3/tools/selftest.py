#!/usr/bin/env python3
"""Source-transcription checks for Quartz-Seq2 oligos."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import quartz_seq2 as Q
from checks import Check
check = Check(); check.section("2018 Supplementary Table S4 transcription")
check("eMDRT0001", Q.RT_PRIMER_1536_FIRST, "TATAGAATTCGCGGCCGCTCGCGATACATCAATCCTATCTGCNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTT")
check("Tagging primer", Q.TAGGING_PRIMER, "TATAGAATTCGCGGCCGCTCGCGATTTTTTTTTTTTTTTTTTTTTTTT")
check("gM_primer", Q.GM_PRIMER, "GTATAGAATTCGCGGCCGCTCGCGAT")
check("rYshapeP5", Q.RYSHAPE_P5, "GATCGGAAGAGCGTCGTGTA")
check("rYshapeP7LT06", Q.RYSHAPE_P7_LT06, "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
check("P5-gMac_hybrid", Q.P5_GMAC, "AATGATACGGCGACCACCGAGATCTACATTGTATAGAATTCGCGGCCGCTCGCGATAC")
check("Read1DropQuartz", Q.READ1_DROPQUARTZ, "ACATTGTATAGAATTCGCGGCCGCTCGCGATAC")
check.report()
