#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sherry
from checks import Check
import nextera as nx
c=Check(); c.section("SHERRY source transcription")
c("ME oligo",nx.ME_RC,"CTGTCTCTTATACACATCT")
c("Adaptor A",nx.S5+nx.ME,"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("Adaptor B",nx.S7+nx.ME,"GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")
c.report()
