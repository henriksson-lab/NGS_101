#!/usr/bin/env python3
"""Transcription checks against the defining original, not the rendered model."""
from pathlib import Path
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from storage import data_root
from checks import Check, Source, have
import chem

check = Check()
bundle = data_root() / 'sources' / HERE.parent.name
source = Source(bundle / 'nmeth.3478_PMC4624326.xml.txt',
                'doi:10.1038/nmeth.3478; ref/MANIFEST.md')
if have(check, source, label='defining Methods original'):
    original = re.sub(r'\s+', '', source.path.read_text())
    check('complete RNA starter including 3-prime block',
          chem.RNA_SEQUENCE + '/3SpC3/' in original, True)
    check('complete multiplex PCR primer', chem.DM_FORWARD.top in original, True)
    reverse = chem.REVERSE_PARTS[0].top + 'BARCODE' + chem.REVERSE_PARTS[-1].top
    check('complete indexed PCR primer with source placeholder', reverse in original, True)
check.report()
