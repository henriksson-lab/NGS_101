#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import fs_pcrfree as M
c=Check()
s=Source(sources_root()/H.parent.name/'manuale7430_e7435.pdf.txt','Original printed vendor revision; see ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('printed source identity','Version2.0_6/25' in text,True)
 for term in ['singlestrandregionsoftheadaptors', 'E7395', 'DTT']:c('source reagent/geometry '+term,term in text,True)
c.report()
