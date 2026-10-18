#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import ramda_seq as m
c=Check()
s=Source(sources_root()/H.parent.name/'s41467-018-02866-0_supp_41467_2018_2866_MOESM6_ESM.xlsx.txt','Published Supplementary Data3, doi:10.1038/s41467-018-02866-0; tools/get_sources.py RamDA-seq')
if have(c,s):
 rows=dict(line.split('\t',1) for line in s.path.read_text().splitlines() if '\t' in line)
 c('first NSR ordering example matches original Data3',m.FIRST_NSR_1,rows['1st_MmNSR_1'])
 c('complementary second NSR matches original Data3',m.SECOND_NSR_1,rows['2nd_MmNSR_1'])
c.report()
