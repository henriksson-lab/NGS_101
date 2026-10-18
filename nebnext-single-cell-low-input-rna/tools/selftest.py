#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import low_input_rna as M
c=Check()
s=Source(sources_root()/H.parent.name/'manuale6420.pdf.txt','Source original manual printed E6420v6.0_01/24; see local ref/MANIFEST.md')
if have(c,s):
 text=''.join(s.path.read_text().split())
 c('manual identity','Version6.0_01/24' in text,True)
 c('printed TSO including riboguanosines',M.TSO in text,True)
 c('printed RT T30V oligo',M.RT_PRIMER in text,True)
 c('printed one-species cDNA PCR primer',M.CDNA_PCR_PRIMER in text,True)
c.report()
