#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H.parents[1]/"lib")]
from checks import Check
from chem import (BAG_DNA_CAPTURE_WRITTEN, BAG_RNA_CAPTURE_WRITTEN,
                  HYBRID_TSO_WRITTEN, CUSTOM_READ1_PRIMER,
                  P5_ADAPTER3, N701_PRIMER)
c=Check(); c.section("hybrid BAG source transcription")
c("DNA capture oligo", BAG_DNA_CAPTURE_WRITTEN,
  "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG NBBBB KKKKKKKGKKKKKKKKNN")
c("RNA capture oligo", BAG_RNA_CAPTURE_WRITTEN,
  "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGG DDDDDNNN T(30)")
c("hybrid template-switch oligo", HYBRID_TSO_WRITTEN,
  "AAGGCTCTTGACATGAAGGATAGrGrGrG")
c("Read1-Adapter3 custom sequencing primer", CUSTOM_READ1_PRIMER,
  "CGCTGTCCAAGGCTCTTGACACGAAGGATAG")
c("P5-Adapter3 final PCR primer", P5_ADAPTER3,
  "AATGATACGGCGACCACCGAGATCTACACCGCTGTCCAAGGCTCTTGACACGAAGGATAG")
c("N701 illustrative final PCR primer", N701_PRIMER,
  "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG")
c.report()
