#!/usr/bin/env python3
"""Check source-specific modified hairpin transcription in an original."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from storage import sources_root
c=Check()
s=Source(sources_root()/"trael-seq__10.1371+journal.pbio.3000886/journal.pbio.3000886_PMC8021198.xml.txt","Published PLOS Methods adapter sequences")
if have(c,s,label="TrAEL published modified hairpins"):
 text=re.sub(r"\s+","",s.path.read_text())
 c("source adapter1 biotin-bearing sequence",P.ADAPTER1.replace("ATTGGCCTGCGC","ATTGGCC[BtndT]GCGC") in text,True)
 c("source adapter2 terminal phosphorothioate",P.ADAPTER2[:-1]+"*"+P.ADAPTER2[-1] in text,True)
c.report()
