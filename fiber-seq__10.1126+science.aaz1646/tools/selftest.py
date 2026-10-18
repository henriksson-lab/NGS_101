#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import fiber_seq as M
c=Check()
s=Source(sources_root()/H.parent.name/'aaz1646_stergachis_sm.pdf.txt','Published original; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split()).lower()
 c('originalSM printedDOI',('10.1126/science.aaz1646').lower() in text,True)
 c('original200Umarking notlivingprotocol',('using200units(2.5μl)ofHia5').lower() in text,True)
 c('sourceAMpurethenmechanicalshearing',('shearedusingaMegaruptor').lower() in text,True)
c.report()
