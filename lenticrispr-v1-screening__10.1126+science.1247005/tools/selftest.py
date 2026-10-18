#!/usr/bin/env python3
"""Original deposited-component transcription; no repeated construction invariants."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import lenticrispr_v1 as P
from checks import Check,Source,have
from storage import legacy_reference,sources_root
from plasmid import read_genbank
c=Check()
path=legacy_reference('lenticrispr-v1-screening__10.1126+science.1247005/ref/plasmids/addgene-49535_lentiCRISPRv1.gb')
original=Source(path,'https://www.addgene.org/49535/ ; external reference manifest records deposited sequence identity')
if have(c,original,label='deposited invariant cassette components'):
    plasmid=read_genbank(str(path))
    u6=next(f for f in plasmid.features if f.label=='U6 promoter')
    scaffold=next(f for f in plasmid.features if f.label=='gRNA scaffold')
    c('reused hU6 terminal component matches original depositedfeature',P.U6,plasmid.sub(u6.start,u6.end)[-len(P.U6):])
    c('reused invariant scaffold matches original depositedfeature',P.SCAFFOLD,plasmid.sub(scaffold.start,scaffold.end))
c.report()
