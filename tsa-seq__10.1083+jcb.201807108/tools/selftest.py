#!/usr/bin/env python3
"""Published TSA library/run numerals without inferring an acquisition mask."""
import re
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import chem as M
from storage import data_root
from checks import Check,Source,have
c=Check()
s=Source(data_root()/'sources'/H.parent.name/'jcb.201807108_PMC6219710.xml.txt','Chen2018 published Methods library preparation')
if have(c,s,label='TSA published library methods'):
 text=s.path.read_text();m=re.search(r'indexed adaptors with (\d+)-nt barcodes',text)
 c('physical adapter-barcode width',int(m.group(1)) if m else None,M.INDEX_ROLE.length_bp)
 m=re.search(r'sequenced for (\d+) cycles',text)
 c('reported sequencing cycles',int(m.group(1)) if m else None,M.REPORTED_CYCLES)
 m=re.search(r'million (\d+)-nt reads',text)
 c('reported output-read length',int(m.group(1)) if m else None,M.REPORTED_READ_NT)
c.report()
