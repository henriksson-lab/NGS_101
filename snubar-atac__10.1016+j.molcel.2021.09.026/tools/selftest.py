#!/usr/bin/env python3
"""Source-transcription boundaries not enforced by molecular construction."""
from pathlib import Path
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from checks import Check, Source, have
from storage import sources_root
import chem

check = Check()
source = Source(sources_root() / HERE.parent.name /
                'j.molcel.2021.09.026_1-s2.0-S1097276521007954-mmc1.pdf.txt',
                'Wang et al. 2021 published Document S1, Figure S1A and Methods S1')
if have(check, source):
    text = source.path.read_text()
    # Figure's I8 notation is intentionally not recognized as N8 by the shared
    # whole-strand text matcher. Check its actual printed family independently.
    match = re.search(r'5-([ACGT]+)-I8-([ACGT]+)-3', text)
    check('Figure S1A supplies an indexed ATAC reverse-primer family', bool(match))
    if match:
        check('ATAC reverse-primer arms and eight-base index agree with original',
              chem.ATAC_REVERSE, match[1] + 'N' * 8 + match[2])
    check('ATAC forward primer is printed in Figure S1A',
          '5-' + chem.ATAC_FORWARD + '-3' in text)
    match = re.search(r'Cell barcodes primer 5-([ACGT]+)-N16-([ACGT]+)-3', text)
    check('Figure S1A supplies cell-barcode primer family', bool(match))
    if match:
        check('GEM cell-barcode primer agrees with printed Figure S1A family',
              chem.protocol_graph().materials['cell_barcode_primer'].state.top(),
              match[1] + 'N' * 16 + match[2])
    match = re.search(r'5-([ACGT]+)-I6-([ACGT]+)-3', text)
    check('Figure S1A supplies RPIX primer family', bool(match))
    if match:
        check('RPIX arm agrees with original primer, preserving figure discrepancy',
              chem.RPIX_PRIMER, match[1] + 'N' * 6 + match[2])
check.report()
