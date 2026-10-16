#!/usr/bin/env python3
"""Whole later comparator transcription; no original2013 oligo identity implied."""
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from storage import data_root
from checks import Check,Source,have
import pro_cap as chem
check=Check()
b=data_root()/'sources'/HERE.parent.name
source=Source(b/'PMC5502525_PMC5502525.xml.txt','Mahat2016 Table1 comparator, not assigned historical2013 reagents')
if have(check,source,label='later2016 original comparator table'):
 text=re.sub(r'\s+','',source.path.read_text())
 check('complete comparator RA3 including inverted-dT block',chem.RA3_2016+'-/InverteddT/' in text,True)
 check('complete comparator RA5 RNA',chem.RA5_2016 in text,True)
 check('complete comparator RTP',chem.RTP_2016 in text,True)
 check('complete comparator RP1',chem.RP1_2016 in text,True)
 check('complete comparator RPI-n family',chem.RPI_2016 in text,True)
check.report()
