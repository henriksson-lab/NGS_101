#!/usr/bin/env python3
"""Source-transcription checks for this protocol's native declarations."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as p
c=Check()
s=Source(sources_root()/H.parent.name/'nbt.2282_PMC3467340.xml.txt','Published defining SMART-seq article, PMC3467340; tools/get_sources.py SMART-seq')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('CDS primer matches original T(30)VN notation',p.DT_SEQUENCE.replace('T'*30,'T(30)') in text,True)
t=Source(sources_root()/H.parent.name/'nbt.2282_41587_2012_BFnbt2282_MOESM12_ESM.xls.txt','Published Supplementary Table1; tools/get_sources.py SMART-seq')
if have(c,t):
 c('original Table1 explicitly states single-end acquisition','sequenced from a single-end of the fragments' in t.path.read_text(),True)

c.report()
