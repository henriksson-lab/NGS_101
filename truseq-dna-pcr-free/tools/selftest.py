#!/usr/bin/env python3
"""Optional original-source transcription checks; construction enforces pairing."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from checks import Check, Source, have
from storage import sources_root
import chem
c=Check()
s=Source(sources_root()/'truseq-dna-pcr-free'/'illumina_1000000002694_22_ilmn-adapter-sequences.pdf.txt',
         'Illumina Adapter Sequences 1000000002694 v22, pp.79–80')
if have(c,s):
    text=''.join(s.path.read_text().split())
    u,i=chem.arms()
    c('v22 CD i5 arm transcription', ''.join(x.top for x in u).replace('J'*8,'[i5]') in text)
    c('v22 CD i7 arm transcription', ''.join(x.top for x in i).replace('I'*8,'[i7]') in text)
chem.protocol_graph()
c.report()
