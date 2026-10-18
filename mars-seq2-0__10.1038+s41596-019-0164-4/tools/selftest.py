#!/usr/bin/env python3
"""Check source-specific transcription against cached original text, not model labels."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check();root=sources_root()/"mars-seq2-0__10.1038+s41596-019-0164-4"
s=Source(root/"s41596-019-0164-4.pdf.txt","published2019protocol reagent setup")
if have(c,s,label="2019native RT2/PCRoligos"):
 text=re.sub(r"\s+","",s.path.read_text())
 c("published5to3RT2",P.RT2 in text,True)
 c("published P5/Read1 enrichment strand",P.P5_PCR in text,True)
 c("published P7/Read2 enrichment strand",P.P7_PCR in text,True)
s=Source(root/"scg_41596_2019_164_MOESM4_ESM.xlsx.txt","originalSupplementTable1 A1")
if have(c,s,label="2019source RT1table"):
 text=re.sub(r"\s+","",s.path.read_text())
 ordered="".join("CTATTCG" if p.name=="well barcode" else p.top for p in P.rt_parts())
 c("source A1 well7/UMI8/T20N RT1",ordered in text,True)
s=Source(root/"scg_41596_2019_164_MOESM5_ESM.xlsx.txt","originalSupplementTable2 ix1")
if have(c,s,label="2019source blockedadapter"):
 text=re.sub(r"\s+","",s.path.read_text())
 c("source physical5P/N5/C3blocked ix1 donor", "/5Phos/"+"".join(p.top for p in P.adapter_parts())+"/3SpC3/" in text,True)
c.report()
