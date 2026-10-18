#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import invadeseq as M
M.protocol_graph()
source=Source(sources_root()/H.parent.name/"s41596-023-00888-7_PMC10790651.xml.txt","Defining2023 published protocol Materials/Methods")
if have(c,source):
    text="".join(source.path.read_text().split())
    for name,value in (("RT16S1100R",M.RT_16S),("EnrichForward",M.ENRICH_FORWARD),("NestedReverse",M.ENRICH_REVERSE)):
        c("published protocol "+name,value in text,True)
c.report()
