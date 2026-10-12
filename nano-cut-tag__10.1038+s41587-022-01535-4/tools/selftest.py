#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import nano_cut_tag as P
from checks import Check
c=Check();c.section("nano-CT oligo transcription")
c("barcode A",P.BARCODE_A,"TATAGCCT")
c("outer and inner connectors",(P.OUTER,P.CONNECTOR),("TCCACGC","GCGATCGAGGACGGC"))
c("custom sequencing primers",(P.CUSTOM_R1,P.CUSTOM_I2),
  ("GCGATCGAGGACGGCAGATGTGTATAAGAGACAG","CTGTCTCTTATACACATCTGCCGTCCTCGATCGC"))
P.workflow();c.report()
