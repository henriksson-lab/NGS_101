#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import gt_seq as M
M.protocol_graph()
src=Source(sources_root()/H.parent.name/"zenodo895966_nmeth.3370_article.pdf.txt","2015 accepted manuscript OnlineMethods, Zenodo895966")
if have(c,src):
    text="".join(src.path.read_text().split())
    c("2015 capture source shorthand",M.CAPTURE_HANDLE+"T30VN" in text,True)
    c("2015 TSO modifications and version",M.TSO_WRITTEN in text,True)
    c("2015 universal PCR source",M.ISPCR+"-3'" in text,True)
c.report()
