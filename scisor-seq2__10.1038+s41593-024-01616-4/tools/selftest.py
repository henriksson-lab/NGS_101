#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as M
c=Check()
s=Source(sources_root()/H.parent.name/'s41593-024-01616-4_PMC11156538.xml.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('nativePartialRead1',(M.PARTIAL_READ1).lower() in text,True)
 c('nativePartialTSO',(M.PARTIAL_TSO).lower() in text,True)
 c('nativeLSK110 guideidentity',('ACDE_9110_v110_revC_10Nov2020').lower() in text,True)
c.report()
