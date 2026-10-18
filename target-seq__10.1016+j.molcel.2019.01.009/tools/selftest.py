#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as M
c=Check()
s=Source(sources_root()/H.parent.name/'j.molcel.2019.01.009_PMC6436961.xml.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('sourceISPCR',(M.ISPCR).lower() in text,True)
 c('physicalTSO RNA/LNA',(M.TSO_BASES[:-3]+'rGrG+G').lower() in text,True)
 c('nativeCS1tag',(M.CS1).lower() in text,True)
 c('nativeCS2tag',(M.CS2).lower() in text,True)
if have(c,s):c('whole source DT IUPACstrand',''.join(p.top for p in M.DT).lower() in text,True)
c.report()
