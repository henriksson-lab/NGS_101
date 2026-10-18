#!/usr/bin/env python3
"""Source transcription; native construction owns geometry invariants."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from checks import Check, Source, have
from storage import sources_root
import chem
c=Check()
s=Source(sources_root()/'yamat-seq__10.1093+nar+gkx005'/'nar_gkx005_PMC5605243.xml.txt',
         'Published PMC5605243 Methods, Y-shaped adapter design')
if have(c,s):
    text=''.join(s.path.read_text().split())
    c('Y3 ordered phosphate/ddC transcription', '5phos/'+chem.Y3_SEQUENCE+'/3ddC' in text)
    for base in 'AGCU':
        c('Y5-'+base+' DNA/RNA transcription', chem.Y5_DNA+'g'+base.lower() in text)
chem.protocol_graph()
c.report()
