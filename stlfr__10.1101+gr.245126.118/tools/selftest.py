#!/usr/bin/env python3
"""Check local stLFR source choices against the original order workbook."""
from pathlib import Path
import os
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'lib'))
import stlfr_chem as model
from checks import Check, Source, have

check = Check()
data = Path(os.environ.get('CHEM_DATA', HERE.parents[1] / '_data'))
source = Source(data / 'sources' / HERE.parent.name / '324392_324392-2.xlsx.txt',
                'Original Supplement_0.xlsx, Cheng2018 doi:10.1038/protex.2018.116; recovered protocols.io attachment7 and bioRxiv324392 supplement2')
if have(check, source, label='Original stLFR order workbook'):
    sheets = {}
    for line in source.path.read_text().splitlines():
        if line.startswith('## sheet '):
            sheet = line.removeprefix('## sheet ')
            sheets[sheet] = []
        else:
            sheets[sheet].append([field.strip() for field in line.split('\t')])
    common = {}
    scope = 'common'
    for row in sheets['Common adapter sequences']:
        if row[0].startswith('For '):
            scope = row[0]
        elif len(row) > 1:
            common[scope, row[0]] = row[1]
    pairs = [('common', 'BeadCommonT', '/52-Bio/' + model.BEAD_COMMON_T),
             ('common', 'Bridge', model.BRIDGE_RNA),
             ('For BGISEQ-500', 'Transposon1T', '/5Phos/' + model.TRANSPOSON_1_T),
             ('For BGISEQ-500', 'TransposonB', '/5Phos/' + model.TRANSPOSON_B),
             ('For BGISEQ-500', 'PCR1', model.PCR1),
             ('For BGISEQ-500', 'PCR2', model.PCR2),
             ("For 3' Branch Ligation", 'BranchT', '/5Phos/' + model.BRANCH_T),
             ("For 3' Branch Ligation", 'BranchB', model.BRANCH_B + "/3'dd/G"),
             ('For Gap Ligation', 'PCR2', model.ILLUMINA_PCR2)]
    check('selected common/branch oligos retain original ends and platform-specific primers',
          [(scope, name, sequence) for scope, name, sequence in pairs],
          [(scope, name, common[scope, name]) for scope, name, _ in pairs])
    for number, prefix, suffix in ((1, 'CTCTGG', 'TCTGCG'), (2, '', ''),
                                    (3, 'CCTTCC', 'TAGCATGGACTATGG')):
        rows = [r for r in sheets[f'Barcode{number} Top'] if len(r) == 4 and r[0] != 'Name']
        check(f'barcode set{number}: complete source members fit declared ten-base architecture',
              len(rows) == 1536 and all(len(r[2]) == 10 and r[1] == '/5Phos/' + prefix + r[2] + suffix for r in rows))
check.report()
