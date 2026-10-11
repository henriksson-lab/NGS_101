#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import chem
from accuracy_protocols import TIME_B_TS301, TIME_CUSTOM_I1, TIME_METHYLATED_A
from checks import Check

# These checks protect transcription from Supplementary Table 5. Primer landing,
# feature spans and workflow continuity are enforced by their constructors.
c = Check()
c.section("TIME-seq Supplementary Table 5 transcription")
c("TS301 indexed adaptor B", TIME_B_TS301,
  "TGGGTGGAGGGTGGAAAGAAGATGTGTATAAGAGACAG")
c("methylated adaptor A, with IDT modification notation", TIME_METHYLATED_A,
  "T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG")
c("custom Index Read 1 primer", TIME_CUSTOM_I1,
  "CCACCCTCCACCCAGATCTGCGTTGATACCACTGCTT")
c.report()
