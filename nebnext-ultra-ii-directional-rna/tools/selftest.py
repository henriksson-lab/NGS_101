#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import directional_rna as M
c=Check()
s=Source(sources_root()/H.parent.name/'manuale7760_e7765.pdf.txt','Original printed vendor revision; see ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('printed source identity','Version7.0_4/25' in text,True)
 for term in ['dUTP', 'StrandSpecificityReagent', 'USER']:c('source reagent/geometry '+term,term in text,True)
c.report()
