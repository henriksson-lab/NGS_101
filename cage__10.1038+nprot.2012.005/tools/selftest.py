#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import cage as C
from checks import Check
c=Check(); c.section("source transcription")
c("RT-N15-EcoP", C.RT_ECOP, "AAGGTCTATCAGCAGNNNNNNNNNNNNNNN")
c("PCR forward", C.FWD, "AATGATACGGCGACCACCGACAGGTTCAGAGTTC")
c("PCR reverse", C.REV, "CAAGCAGAAGACGGCATACGA")
c("custom sequencing primer", C.SEQUENCING, "CGGCGACCACCGACAGGTTCAGAGTTCTACAG")
c.report()
