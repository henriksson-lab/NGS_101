#!/usr/bin/env python3
"""Whole-source transcription checks; construction validates primer landing."""
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from storage import data_root
from checks import Check,Source,have
import pro_seq
check=Check()
bundle=data_root()/'sources'/HERE.parent.name
source=Source(bundle/'594009_926d2d49ebf94caa8e6cf707c2f4846a.pdf.txt','Original2013 SOM author-hosted copy; ref/MANIFEST.md')
if have(check,source,label='original2013 SOM'):
    text=re.sub(r'\s+','',source.path.read_text()).replace('′',"'").replace('’',"'")
    rna=lambda dna: ''.join('r'+base for base in dna.replace('T','U'))
    check('complete reverse3 RNA adapter including both modifications',
          "5'p-"+rna(pro_seq.VRA3)+"-/3'InvdT/" in text,True)
    check('complete original TruSeq reverse5 RNA option',rna(pro_seq.VRA5) in text,True)
    check('complete original RP1 primer',pro_seq.RP1 in text,True)
source=Source(bundle/'PMC5502525_PMC5502525.xml.txt','Mahat2016 Table1 comparator, not original2013 kit identity')
if have(check,source,label='later explicit RPI-n comparator'):
    text=re.sub(r'\s+','',source.path.read_text())
    check('complete later2016 RPI-n ordering family',pro_seq.RPI in text,True)
check.report()
