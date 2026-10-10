#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tgirt_seq as M
from checks import Check
c=Check(); c.section("TGIRT-seq source transcription")
c("R2 RNA as ordered",M.R2_RNA_WRITTEN,"rArArGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrUrCrArC/3SpC3/")
c("R2R DNA",M.R2R_DNA,"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTN")
c("R1R DNA as ordered",M.R1R_WRITTEN,"/5Phos/GATCGTCGGACTGTAGAACTCTGAACGTGTAG/3SpC3/")
c.report()
