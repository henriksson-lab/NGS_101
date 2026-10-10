#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import snmc2t_seq as P
from checks import Check
c=Check();c.section("snmC2T modality invariant");c("one mixed final library",P.FINAL_LIBRARY.name,"snmC2T-seq mixed DNA/cDNA library");c.report()
