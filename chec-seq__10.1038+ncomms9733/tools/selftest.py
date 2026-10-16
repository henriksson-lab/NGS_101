#!/usr/bin/env python3
"""Check the illustrated Orsi variant against its actual source oligos."""
import os
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from checks import Check, Source, have
from chemdraw import revcomp
import chem

check = Check()
data = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data'))
source = Source(data / 'sources' / HERE.parent.name / 'PMC4410783_PMC4410783.xml.txt',
                'Original Orsi et al. Basic Protocol3, doi:10.1002/0471142727.mb2131s110; PMC4410783, adapter/PCR oligos812–832')
if have(check, source, label='Original Orsi adapter and PCR oligos'):
    text = source.path.read_text()
    universal = re.search(r'TruSeq Universal Adapter:\s*5[’\x27]-([ACGTN*]+)', text)[1]
    indexed = re.search(r'TruSeq Indexed Adapter:\s*5[’\x27]-\[P\]([ACGTN*]+)', text)[1]
    p5 = re.search(r'P5 PCR Primer:\s*5[’\x27]\s*([ACGT]+)', text)[1]
    p7 = re.search(r'P7 PCR Primer:\s*5[’\x27]\s*([ACGT]+)', text)[1]
    check('complete adapter/PCR bases agree with original source',
          (chem.CHEC_UNIVERSAL_ADAPTER, chem.CHEC_INDEXED_ADAPTER, chem.CHEC_P5_PCR, chem.CHEC_P7_PCR),
          (universal.replace('*', ''), indexed.replace('*', ''), p5, p7))
    check('displayed adapters retain original terminal phosphorothioate notation',
          universal in chem.oligos()[0] and indexed in chem.oligos()[1])
    library, primers = chem.chec_library()
    check('left endpoint contains the entire source Universal Adapter once',
          library.segments[0].top, universal.replace('*', ''))
    check('endpoint retains the complete indexed adapter including source six-position index',
          ''.join(s.top for s in library.segments[-3:]), indexed.replace('*', ''))
    r1 = chem.sp.locate(library, chem.sp.TRUSEQ['R1'])
    check('Read1 starts immediately after the whole source Universal Adapter',
          r1.end, len(universal.replace('*', '')))
    check('both source enrichment primers have full sites on the amplified endpoint',
          library.top().startswith(p5) and revcomp(p7) in library.top())
check.report()
