#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import cite_seq as P
from checks import Check,Source,have
from storage import sources_root
from multimodal_page import render
c=Check();render(P)
source=Source(sources_root()/H.parent.name/"nmeth.4380_41592_2017_BFnmeth4380_MOESM3_ESM.pdf.txt","published author protocol version2017-07-11")
if have(c,source,label="author protocol modified oligos"):
 text="".join(source.path.read_text().split())
 for name,seq in (("v2_BC1 5′amine/C12",P.V2_BC1),("PE1.0 P5 phosphorothioates",P.PE1_P5),("RPI1 index primer",P.RPI1)):
  c("published "+name,seq in text,True)
c.report()
