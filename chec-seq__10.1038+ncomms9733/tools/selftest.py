#!/usr/bin/env python3
"""Check the illustrated Orsi variant against its actual source oligos."""
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from storage import data_root
from checks import Check, Source, have
from chemdraw import revcomp
import chem

check = Check()
data = data_root()
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
check.report()
