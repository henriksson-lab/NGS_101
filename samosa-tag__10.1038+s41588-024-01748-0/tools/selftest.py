#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import samosa_tag as M
c=Check()
s=Source(sources_root()/H.parent.name/'s41588-024-01748-0_supp_41588_2024_1748_MOESM7_ESM.xlsx.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('complete modifiedbc001 hairpin',('/5Phos/'+M.ADAPTER).lower() in text,True)
c.report()
