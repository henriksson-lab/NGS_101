#!/usr/bin/env python3
"""Transcribe original Table1; constructor functions enforce library assembly."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import moffat_readout as P
from checks import Check,Source,have
from storage import sources_root
c=Check();original=Source(sources_root()/H.parent.name/'moffat_lab_crispr_screen_protocol_v2.pdf.txt','Author protocol, printedREV.20190723 Table1')
if have(c,original,label='original Table1'):
 text=original.path.read_text().split('Table 1. PCR primers',1)[1].split('Table 2.',1)[0]
 for label,value in [('LCV2 forward primer',P.LCV2_PCR1_F),('LCV2 reverse primer',P.LCV2_PCR1_R),('D501-F',P.D501_F),('D701-R',P.D701_R)]:
  hit=re.search(re.escape(label)+r'\s+([ACGT \n]+)',text)
  c(label+' transcription',value,re.sub(r'\s+','',hit.group(1)) if hit else None)
c.report()
