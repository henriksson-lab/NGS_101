#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as M
c=Check()
s=Source(sources_root()/H.parent.name/'illumina_15031942_07_nextera-xt-product-documentation.pdf.txt','Original printed vendor revision; see ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('printed source identity','15031942v07' in text,True)
 for term in ['0.1NNaOH', 'LNW1', 'LNS1']:c('source reagent/geometry '+term,term in text,True)
c.report()
