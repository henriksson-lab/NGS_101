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
s=Source(sources_root()/H.parent.name/'s41587-020-0497-0_41587_2020_497_MOESM3_ESM.xlsx.txt','Published SMART-seq3 Supplementary Table1; tools/get_sources.py SMART-seq3')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('native biotinylated dT primer from published Table1','/5BiosG/'+p.DT_SEQUENCE in text,True)
 c('native biotinylated N8 RNA-tail TSO from published Table1','/5Biosg/'+p.TSO_HANDLE+'N'*8+'rGrGrG' in text,True)
 c('published asymmetric forward primer',p.FWD in text,True)
t=Source(sources_root()/H.parent.name/'protocolsio_smart-seq3-protocol-36wgq5rjxgk5.pdf.txt','Author protocol V.3, printed bcq4ivyw; tools/get_sources.py SMART-seq3')
if have(c,t):
 text=''.join(t.path.read_text().split())
 c('native protected forward primer',p.FWD[:-3]+p.FWD[-3]+'*'+p.FWD[-2]+'*'+p.FWD[-1] in text,True)
 c('native protected reverse primer',p.REV[:-3]+p.REV[-3]+'*'+p.REV[-2]+'*'+p.REV[-1] in text,True)

c.report()
