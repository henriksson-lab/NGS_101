#!/usr/bin/env python3
"""Check complete loaded oligos against the cited original, including phosphate."""
import os
from pathlib import Path
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from checks import Check, Source, have
import cut_tag

check = Check()
bundle = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data')) / 'sources' / HERE.parent.name
source = Source(bundle / 'PMC4248319_PMC4248319.xml.txt',
                'Picelli2014 doi:10.1101/gr.177881.114; ref/MANIFEST.md')
if have(check, source, label='Picelli transposome assembly original'):
    original = re.sub(r'\s+', '', source.path.read_text())
    check('complete source ME-A', '5′-' + cut_tag.ME_A + '-3′' in original, True)
    check('complete source ME-B', '5′-' + cut_tag.ME_B + '-3′' in original, True)
    check('complete source MErev with phosphate',
          '5′-[phos]' + cut_tag.ME_REVERSE + '-3′' in original, True)
check.report()
