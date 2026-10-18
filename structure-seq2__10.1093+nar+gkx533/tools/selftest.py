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
s=Source(sources_root()/H.parent.name/'nar_gkx533_PMC5737731.xml.txt','Published Structure-seq2 article PMC5737731; tools/get_sources.py Structure-seq2')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('source phosphorylated/C3 blocked hairpin','p'+p.SS2_HAIRPIN_DONOR.top()+'-3′-Spacer' in text,True)
 c('published forward primer',p.SS2_FORWARD_SEQUENCE in text,True)
 c('published reverse fixed arm',p.SS2_REVERSE_SEQUENCE[-33:] in text,True)
# Figure-only RT primer is externally transcribed from rendered Supplementary FigureS1.
c('FigureS1 RT primer transcription',p.SS2_RT_SEQUENCE,'CAGACGTGTGCTCTTCCGATCNNNNNN')

c.report()
