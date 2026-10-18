#!/usr/bin/env python3
"""Pin the defining2012 times without silently using the2013 oxidation variant."""
import re
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import chem as M
from storage import data_root
from checks import Check,Source,have
c=Check()
s=Source(data_root()/'sources'/H.parent.name/'j.cell.2012.04.027_PMC3589129.xml.txt','Yu2012 defining ExperimentalProcedures')
if have(c,s,label='TAB defining2012 reaction methods'):
 text=s.path.read_text();text=text[text.index('Glucosylation and Oxidation of Genomic DNA'):]
 times=re.findall(r'reaction was incubated at 37\s*°C for ([\d.]+) h',text)
 c('defining glucosylation/oxidation minutes',[int(float(t)*60) for t in times[:2]],[M.GLUCOSYLATION_MIN,M.OXIDATION_2012_MIN])
c.report()
