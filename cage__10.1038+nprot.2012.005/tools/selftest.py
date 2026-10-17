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
c("5-prime SOL N6 AGA", C.FIVE_UPPER_N6,
  "CCACCGACAGGTTCAGAGTTCTACAGAGACAGCAGNNNNNN")
c("5-prime SOL GN5 AGA", C.FIVE_UPPER_GN5,
  "CCACCGACAGGTTCAGAGTTCTACAGAGACAGCAGGNNNNN")
c("5-prime SOL lower AGA", C.FIVE_LOWER,
  "CTGCTGTCTCTGTAGAACTCTGAACCTGTCGGTGG")
c("3-prime upper linker", C.THREE_UPPER, "NNTCGTATGCCGTCTTCTGCTTG")
c("3-prime lower linker", C.THREE_LOWER, "CAAGCAGAAGACGGCATACGA")
c.section("model construction")
C.protocol_graph()
c.report()
