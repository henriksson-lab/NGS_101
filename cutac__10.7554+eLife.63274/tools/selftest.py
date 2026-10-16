#!/usr/bin/env python3
"""Whole original-source oligos, including source modification and family arms."""
import json
import os
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from checks import Check, Source, have
import cutac
check=Check()
base=Path(os.environ.get('CHEM_DATA',HERE.parents[1]/'_data'))/'sources'
source=Source(base/HERE.parent.name/'protocolsio_29425_public.json','doi:10.17504/protocols.io.8yrhxv6; ref/MANIFEST.md')
def source_strings(value):
    if isinstance(value,dict):
        for child in value.values(): yield from source_strings(child)
    elif isinstance(value,list):
        for child in value: yield from source_strings(child)
    elif isinstance(value,str):
        yield re.sub(r'\s+','',value)
        if value.startswith(('{','[')):
            try: nested=json.loads(value)
            except ValueError: return
            yield from source_strings(nested)
if have(check,source,label='original loading protocol API table'):
    strings=set(source_strings(json.loads(source.path.read_text())))
    check('complete source ME-A table cell',cutac.ME_A in strings,True)
    check('complete source ME-B table cell',cutac.ME_B in strings,True)
    check('complete source phosphorylated MErev table cell','[Phos]'+cutac.ME_REVERSE in strings,True)
source=Source(base/'doi_10.1038+nature14590'/'nature14590_41586_2015_BFnature14590_MOESM36_ESM.xlsx.txt','doi:10.1038/nature14590 original MOESM36')
if have(check,source,label='cited original custom-primer workbook'):
    rows=[line.split('\t') for line in source.path.read_text().splitlines()]
    for prefix,pool in [('v2_Ad1.',cutac.P5_POOL),('v2_Ad2.',cutac.P7_POOL)]:
        sequences=[row[1] for row in rows if len(row)>1 and row[0].startswith(prefix)]
        before,after=pool.top().split('N'*8)
        check(prefix+' complete printed family architecture',bool(sequences) and all(len(seq)==len(before)+8+len(after) and seq.startswith(before) and seq.endswith(after) for seq in sequences),True)
check.report()
