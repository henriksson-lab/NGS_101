#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import music as M
from checks import Check
c=Check();c("RNA-linker fixed prefix",M.RNA_LINKER[:12],"CGAGGAGCGCTT")
c("barcode-1 top member 1",M.BC1_TOP,"CTCCTCGGCAACCATTACCTT")
c("universal PCR adapter",M.UNIVERSAL,"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c.report()
