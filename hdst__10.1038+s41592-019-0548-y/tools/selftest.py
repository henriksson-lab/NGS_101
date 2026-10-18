#!/usr/bin/env python3
"""Source transcription against cached originals; source absence skips these checks."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
s=Source(sources_root()/H.parent.name/"s41592-019-0548-y_PMC6765407.xml.txt","Published HDST Bead production121")
if have(c,s,label="originalHDST assembled-primer components"):
 text=s.path.read_text();flat=re.sub(r"\s+","",text)
 c("source-written precursor skeleton",P.PRECURSOR_WRITTEN in flat,True)
 c("both source completed bridges",P.B12 in text and P.B34 in text,True)
 first=int(re.search(r"a (\d+) bp .Spatial barcode pool.",text)[1])
 second,third=map(int,re.search(r"adding (\d+) bp and (\d+) bp pools",text).groups())
 umi=int(re.search(r"a (\d+) bp unique molecular identifier",text)[1])
 parts=P.identity_parts()
 c("physical source barcode/UMIwidths, distinct from preprocessing",tuple(len(parts[i].top) for i in (0,2,4,5)),(first,second,third,umi))
c.report()
