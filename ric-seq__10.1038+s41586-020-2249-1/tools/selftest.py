#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'s41586-020-2249-1_41586_2020_2249_MOESM12_ESM.xlsx.txt','published2020 Table9')
if have(c,src,label='published2020 primer family'):
 text=re.sub(r'\s+','',src.path.read_text())
 c('PE1 source transcription',chem.PE1 in text,True)
 for i,index in enumerate(chem.INDEXES,1):c('index primer'+str(i)+' source transcription',chem.reverse_primer(index).top() in text,True)
c('PCR template retains actual postUSER ss population',g.steps['library_pcr'].geometry.template is g.materials['user_surviving'].state,True)
c.report()
