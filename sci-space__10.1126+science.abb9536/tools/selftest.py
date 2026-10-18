#!/usr/bin/env python3
"""Source transcription against cached originals; source absence skips these checks."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
root=sources_root()/H.parent.name
s=Source(root/"protocolsio_creating-sci-space-grids-for-spatial-barcoding-5qpvoy58bg4o.pdf.txt","Cited bm64k9gw,13Feb2021,ID42940,step7.1")
if have(c,s,label="cited2021gridhash source"):
 text=re.sub(r"\s+","",s.path.read_text())
 c("cited DOI identity preserved despite normalized cache name","10.17504/protocols.io.bm64k9gw" in text and "ProtocolIntegerID:42940" in text and "LastModified:February13,2021" in text,True)
 c("source hash scaffold incl Banchor/A32",P.HASH_PREFIX+"[10bp-barcodesequence]-B"+"A"*32 in text,True)
s=Source(root/"caltech_abb9536-Srivatsan-SM.pdf.txt","Published2021SupplementMethodspp4–5")
if have(c,s,label="study-native identifier widths"):
 text=re.sub(r"\s+","",s.path.read_text())
 umi,rt=map(int,re.search(r"Read1:(\d+)cycles\((\d+)base-pairUMIand(\d+)base-pairRTbarcode",text).groups()[1:])
 c("UMI/RTwidths from published native read declaration",(len(P.rt_end()[1].top),len(P.rt_end()[2].top)),(umi,rt))
 idx=[int(re.search(r"Index"+str(i)+r":\d+cycles\((\d+)base-pairPCRindex"+str(i),text)[1]) for i in (1,2)]
 c("native PCRindex widths explicitly stated separately fromcycles",[len(s.top) for s in P.FINAL_LIBRARY if s.feature in (P.IDX5,P.IDX7)],idx)
c.report()
