#!/usr/bin/env python3
"""Optional transcription against the selected published original."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
s=Source(sources_root()/H.parent.name/'s41596-018-0058-x_PMC6279253.xml.txt','Selected published original source; see native migration bounded record')
if have(c,s,label='Original source transcription'):
 compact=re.sub(r'\s+','',s.path.read_text())
 for name,seq in [('GP6',P.GP6.top()),('GP12N9diversity',P.GP12.top()),('PE1Sa',P.PE1SA.top()),('PE1Sb',P.PE1SB.top()),('R2PCR',P.R2_PCR.top())]:c(name+' sourceprimertranscription',seq in compact,True)
c.report()
