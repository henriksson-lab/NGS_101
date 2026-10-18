#!/usr/bin/env python3
"""Transcription checks for the explicitly inherited2015G&T primer reference."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as p
c=Check()
s=Source(sources_root()/'g-t-seq__10.1038+nmeth.3370'/'zenodo895966_nmeth.3370_article.pdf.txt','Published G&T-seq2015 article, doi:10.1038/nmeth.3370, author-original from Zenodo895966; tools/get_sources.py G&T-seq')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('inherited modified TSO matches cited2015 original',p.TSO_WRITTEN in text,True)
 c('inherited capture primer matches originalT30VN shorthand',p.SMART_HANDLE+p.CAPTURE_LINKER+'T30VN' in text,True)
 c('inherited ISPCR matches original written primer',"5'-"+p.ISPCR+"-3'" in text,True)
c.report()
