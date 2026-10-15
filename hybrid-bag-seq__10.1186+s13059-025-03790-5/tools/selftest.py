#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H.parents[1]/"lib")]
from checks import Check
from chem import BAG_DNA_CAPTURE_WRITTEN, BAG_RNA_CAPTURE_WRITTEN, CUSTOM_READ1_PRIMER
c=Check(); c.section("hybrid BAG source transcription")
c("DNA capture oligo", BAG_DNA_CAPTURE_WRITTEN,
  "/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN")
c("RNA capture oligo", BAG_RNA_CAPTURE_WRITTEN,
  "/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
c("Read1-Adapter3 custom sequencing primer", CUSTOM_READ1_PRIMER,
  "CGCTGTCCAAGGCTCTTGACACGAAGGATAG")
c.report()
