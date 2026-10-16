#!/usr/bin/env python3
"""Whole printed first-round tails; primer/circle geometry is constructed locally."""
from pathlib import Path
import os
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from checks import Check,Source,have
import chem
check=Check()
b=Path(os.environ.get('CHEM_DATA',HERE.parents[1]/'_data'))/'sources'/HERE.parent.name
s=Source(b/'4C-seq-from-beginning-to-end_-A-detailed-protocol-for-sample-preparation-and-data-analysis.pdf.txt','Original Krijger Sections3.4/Table1; ref/MANIFEST.md')
if have(check,s,label='original 4C source'):
 text=re.sub(r'\s+','',s.path.read_text())
 check('complete printed reading tail',chem.READING_TAIL in text,True)
 check('complete printed non-reading tail',chem.NONREADING_TAIL in text,True)
check.report()
