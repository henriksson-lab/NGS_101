#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as M
from checks import Check,Source,have
from storage import sources_root
c=Check();M.protocol_graph()
source=Source(sources_root()/H.parent.name/"s41467-019-11049-4_PMC6635368.xml.txt","Published RAGE Methods Targetedcapture and demultiplexing")
if have(c,source):
    text="".join(source.path.read_text().split())
    for name,value in (("sourceR1",M.R1_PRIMER),("sourceTSO",M.TSO_PRIMER)):
        c(name+" transcription",value in text,True)
c.report()
