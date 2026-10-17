#!/usr/bin/env python3
"""Check complete loaded oligos against the cited original, including phosphate."""
from pathlib import Path
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from storage import data_root
from checks import Check, Source, have
import cut_tag

check = Check()
bundle = data_root() / 'sources' / HERE.parent.name
source = Source(bundle / 'PMC4248319_PMC4248319.xml.txt',
                'Picelli2014 doi:10.1101/gr.177881.114; ref/MANIFEST.md')
if have(check, source, label='Picelli transposome assembly original'):
    original = re.sub(r'\s+', '', source.path.read_text())
    check('complete source ME-A', '5′-' + cut_tag.ME_A + '-3′' in original, True)
    check('complete source ME-B', '5′-' + cut_tag.ME_B + '-3′' in original, True)
    check('complete source MErev with phosphate',
          '5′-[phos]' + cut_tag.ME_REVERSE + '-3′' in original, True)
source = Source(bundle / '41586_2015_BFnature14590_MOESM36_ESM.xlsx.txt',
                'doi:10.1038/nature14590 original MOESM36; catalogue/extra_sources.tsv')
if have(check, source, label='cited original custom-primer workbook'):
    rows = [line.split('\t') for line in source.path.read_text().splitlines()]
    sequences = [row[1] for row in rows if len(row) > 1 and row[0].startswith('v2_Ad2.')]
    before, after = cut_tag.P7_POOL.top().split('N' * 8)
    check('v2_Ad2. complete printed family architecture',
          bool(sequences) and all(len(seq) == len(before) + 8 + len(after)
                                  and seq.startswith(before) and seq.endswith(after)
                                  for seq in sequences), True)
check.report()
