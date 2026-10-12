#!/usr/bin/env python3
"""Source-transcription checks for the published TRIP oligos."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

from checks import Check
import trip as T

c = Check()
c.section("Akhtar et al. 2014 reagent-list transcription")
c("PB-barcode-long-7", T.BARCODE_LONG,
  "GTGACACCTGCAGGATCA" + "N" * 16 + "CTCGAGTTGTGGCCGGCCCTTGTGACTG")
c("indexed forward primer example", T.PB_CDNA_F,
  "ACACTCTTTCCCTACACGACGCTCTTCCGATCTACGAGTGCGTGTCACAAGGGCCGGCCACAA")
c("inverse-PCR primer", T.INV_PCR_F_TRUSEQ_PE2,
  "GATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTACGTCACAATATGATTATCTTTCTAG")
c.report()
