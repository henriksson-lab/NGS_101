#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from storage import sources_root
c=Check();P.protocol_graph()
src=Source(sources_root()/H.parent.name/"nbt.3129_PMC4374170.xml.txt","published original Methods P1/P2/biotinP3")
if have(c,src,label="DR-seq original primer transcription"):
 text="".join(src.path.read_text().split())
 for name,seq in (("P1",P.P1),("P2",P.P2),("biotin-P3",P.P3)):
  c("original "+name,seq in text,True)
c.report()
