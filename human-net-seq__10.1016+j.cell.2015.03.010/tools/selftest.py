#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import human_net_seq as N
from checks import Check
c=Check(); c.section("source transcription")
c("barcode linker", N.LINKER, "NNNNNNCTGTAGGCACCATCAAT")
c("reverse PCR primer", N.REVERSE_PCR, "CAAGCAGAAGACGGCATACGA")
c("custom sequencing primer", N.SEQUENCING, "TCCGACGATCATTGATGGTGCCTACAG")
c.report()
