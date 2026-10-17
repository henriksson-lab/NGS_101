#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
import chem
from checks import Check, Source, have

c = Check()
paper = Source(data_root() / 'sources/sc-end5-seq__10.1101+2021.04.04.438388'
               / "2021.04.04.438388_v2.full.pdf.txt",
               "bioRxiv doi:10.1101/2021.04.04.438388 v2, Methods pp. 31–32")
if have(c, paper, label="source transcription of sc-end5 RT primers"):
    methods = paper.path.read_text(encoding="utf-8")
    c("random RT primer appears in original Methods", chem.RANDOM_RT_PRIMER in methods)
    c("dT primer handle appears in original Methods", chem.RT_HANDLE + "–" in methods)
    c("dT(30)-VN shorthand appears in original Methods", "T(30)–VN" in methods)
    c("13-nt TS-oligo analysis sequence appears in original Methods",
      "(" + chem.TS_OLIGO_ANALYSIS_3PRIME + ")" in methods)
chem.protocol_graph()
c.report()
