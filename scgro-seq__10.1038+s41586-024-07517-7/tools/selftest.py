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
s=Source(sources_root()/H.parent.name/'s41586-024-07517-7_supp_41586_2024_7517_MOESM1_ESM.pdf.txt','Selected published original source; see native migration bounded record')
if have(c,s,label='Original source transcription'):
 compact=re.sub(r'\s+','',s.path.read_text())
 c('design3 physicalazideT/SpC3 withsourceN5/J10',('[azideT]ACAGG'+'N'*5+'J'*10+P.ARM+'[SpC3]') in compact,True)
c.report()
