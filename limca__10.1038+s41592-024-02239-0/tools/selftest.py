#!/usr/bin/env python3
"""Optional transcription checks against original pex2423v1 attachment."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
s=Source(sources_root()/H.parent.name/'3b3ac9309196696067b340eb.xlsx.txt','LiMCA pex2423v1 original attached LiMCAOligos workbook')
if have(c,s,label='LiMCA original oligos'):
 text=s.path.read_text()
 for name,seq in [('META20_1',P.META20_1),('META40_1',P.META40_1),('META40_2',P.META40_2),('ME19',P.ME19),('UDP5example',P.UDP5),('UDP7example',P.UDP7),('sourceTSO',P.TSO_DNA+'rGrG+G')]:c(name+' matches original workbook',seq in text,True)
c.report()
