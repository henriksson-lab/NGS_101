#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import low_bias_small_rna as M
c=Check()
s=Source(sources_root()/H.parent.name/'manuale3420.pdf.txt','Source original manual printed E3420v2.0_5/26; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('manual identity','Version2.0_5/26' in text,True)
 c('source Read1 trim motif',M.MANUAL_TRIM['Read 1'] in text,True)
 c('source Read2 trim motif',M.MANUAL_TRIM['Read 2'] in text,True)
 c('source simultaneous adapter modification creates RT primer','3´adaptorisenzymaticallymodifiedtocreateaprimerforreverse' in text,True)
c.report()
