#!/usr/bin/env python3
"""Source transcription: construction/pairing are enforced by native model."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();src=Source(sources_root()/'scsplat__10.1101+2021.10.14.464375'/'s41598-022-09798-2_supp_41598_2022_9798_MOESM1_ESM.pdf.txt','Published SupplementaryTableS1')
if have(c,src,label='Published scSPLAT oligo table'):
 t=src.path.read_text();c('Blocked representative SSO_H1',re.search(r'SSO_H_1\s+5AmMC6/([ACGTH]+)',t)[1],chem.SSO_TAIL+chem.SSO_H_1_BARCODE+chem.H8)
 c('ss1a arm plus random hexamer',re.search(r'ss1a\s+5[’\']AmMC6/([ACGTN]+)/3[’\']AmMo',t)[1],chem.SS1A)
 c('ss1b phosphate and amine block',re.search(r'ss1b\s+5[’\']Phos/([ACGT]+)/3[’\']AmMo',t)[1],chem.SS1B)
c.report()
