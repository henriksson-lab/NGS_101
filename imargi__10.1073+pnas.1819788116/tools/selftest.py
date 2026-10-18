#!/usr/bin/env python3
"""Optional source transcription, including physical linker modifications."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
s=Source(sources_root()/'doi_10.1038+s41596-019-0229-4'/'s41596-019-0229-4_PMC7314528.xml.txt','Wu2019 published detailed protocol original JATS')
if have(c,s,label='iMARGI modified original oligos'):
 text=s.path.read_text();compact=re.sub(r'\s+','',text)
 c('TOP retains original PS and internalbiotin',P.TOP.replace('TCCGATCT','TCC*GA/iBiodT/CT') in compact,True)
 c('BOTTOM source bases',P.BOTTOM in compact,True)
 c('cut oligo physical inverted terminaldT',P.CUT+'/3InvdT/' in compact,True)
 c('universal primer source terminalPS',P.UNIVERSAL[:-1]+'*'+P.UNIVERSAL[-1] in compact,True)
c.report()
