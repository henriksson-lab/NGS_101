#!/usr/bin/env python3
"""Check source-specific transcription against cached original text, not model labels."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import native as P
c=Check();s=Source(sources_root()/"cel-seq__10.1016+j.celrep.2012.08.003/j.celrep.2012.08.003_mmc4.pdf.txt","original2012CEL-seq detailed protocol RT primer1")
if have(c,s,label="2012CEL-seq primer transcription"):
 text=re.sub(r"\s+","",s.path.read_text())
 # #1 is a source-listed concrete member, independently fixing the model's placeholder window.
 ordered="".join("CATCACGC" if p.name=="cell barcode" else p.top for p in P.rt_parts())
 c("original primer1 including spacer/promoter/RA5/barcode/T24V",ordered in text,True)
c.report()
