#!/usr/bin/env python3
"""Source inventory checks; molecular construction validates topology during render."""
from pathlib import Path
import build_page
from checks import Check, Source, have
from storage import sources_root

build_page.render(build_page.M)
check = Check()
source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
                'IsoSeq_v2_primers_12.fasta',
                'PacBio official IsoSeq_v2_primers_12.fasta, linked from the source note')
if have(check, source, label='official twelve-barcode PCR primer inventory'):
    rows = {}
    name = None
    for line in source.path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].strip()
            rows[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('sequence precedes its FASTA identity')
            rows[name] += line.strip()
    expected = {f'IsoSeqX_bc{i:02d}_5p': build_page.M.barcoded_primer(i)
                for i in range(1, len(build_page.M.ISOSEQ_BARCODES) + 1)}
    expected['IsoSeqX_3p'] = build_page.M.ISOSEQ_3P
    check('all displayed PCR primers agree with complete official FASTA records', expected, rows)
check.report()
