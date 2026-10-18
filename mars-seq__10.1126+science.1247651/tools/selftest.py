#!/usr/bin/env python3
"""Check source-specific transcription against cached original text, not model labels."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import mars_graph as P
c=Check();s=Source(sources_root()/"mars-seq__10.1126+science.1247651/scg_jaitin-sm.pdf.txt","original2014Science SupplementTableS7/S8 extract")
if have(c,s,label="original2014MARS oligos"):
 text=re.sub(r"\s+","",s.path.read_text()).replace("PCRforward","").replace("PCRreverse","")
 c("source RT1 including six X barcode and four N UMI", "".join(p.top for p in P.rt_parts()) in text.replace("X","N"),True)
 c("printed2014RT2 retained",P.RT2 in text,True)
 c("source P5/Read1 primer assembled across role-label line",P.P5_PCR in text,True)
 c("source P7/Read2 primer assembled across role-label line",P.P7_PCR in text,True)
 c("physical5P/C3blocked ix1 donor", "/5Phos/"+"".join(p.top for p in P.adapter_parts())+"/3SpC3/" in text,True)
c.report()
