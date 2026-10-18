#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as M
c=Check()
s=Source(sources_root()/H.parent.name/'s41596-024-01021-y_PMC11537836.xml.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('Smart3 nativebiotin/riboTSO',('biotin-'+''.join(p.top for p in M.TSO[:-1])+'rGrGrG').lower() in text,True)
 c('native forwardPS',(M.FWD[:-3]+'A*T*G').lower() in text,True)
 c('matrix oligoT30V',('ACGAGCATCAGCAGCATACGAT30V').lower() in text,True)
c.report()
