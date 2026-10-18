#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import emseq as M
c=Check()
s=Source(sources_root()/H.parent.name/'manuale8015.pdf.txt','Source original manual printed E8015v2.0_7/26; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('manual identity','Version2.0_7/26' in text,True)
 c('source Read1 trim motif',M.TRIM_READ1 in text,True)
 c('source Read2 trim motif',M.TRIM_READ2 in text,True)
 c('selected low-input branch adds carrier','Add1µlofthe' in text and 'NEBNextCarrierDNAto27µl' in text,True)
c.report()
