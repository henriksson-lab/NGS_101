#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as M
c=Check()
s=Source(sources_root()/H.parent.name/'s41586-019-1367-0_supp_NIHMS1531200-supplement-sup_tab_5.xlsx.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('sourceP5generic',(M.P5).lower() in text,True)
 c('sourceRPIprimer',(__import__('illumina').P7+'XXXXXXXX'+M.RPI).lower() in text,True)
 c('selectedCALRreverse',(M.CALR).lower() in text,True)
c.report()
