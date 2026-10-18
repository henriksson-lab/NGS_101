#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'s41592-019-0570-0_supp_NIHMS1537259-supplement-2.pdf.txt','published Supplementary Protocol C')
if have(c,src,label='published typical sequencing run'):
    text=''.join(src.path.read_text().split()).lower()
    c('selected source SR read length','sr'+str(chem.READ_LENGTHS['Read 1'])+'bpsequencing' in text,True)
c.report()
