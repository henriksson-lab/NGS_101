#!/usr/bin/env python3
"""Deliberately refresh the licensed restriction-data snapshot; no runtime dependency."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
REVISION = 'd7e4b8b19399668b09442a5b35765d9186b5f665'
SOURCE_SHA256 = '2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b'
BASE = f'https://raw.githubusercontent.com/biopython/biopython/{REVISION}/'
SOURCE = BASE + 'Bio/Restriction/Restriction_Dictionary.py'
LICENSE = BASE + 'LICENSE.rst'


def catalogue(source, license_text):
    """Read literal data assignments only; never execute downloaded Python."""
    if hashlib.sha256(source.encode()).hexdigest() != SOURCE_SHA256:
        raise ValueError('Dictionary does not match the pinned release')
    rows = {}
    for node in ast.parse(source).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if (isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name)
                and target.value.id == 'rest_dict'):
            name, value = ast.literal_eval(target.slice), ast.literal_eval(node.value)
            motifs = value['site'].split('|')
            if not all(re.fullmatch('[ACGTRYSWKMBDHVN]+', s) for s in motifs):
                raise ValueError('Unsupported recognition alphabet: ' + name)
            pairs = []
            for a, b in (('fst5','fst3'), ('scd5','scd3')):
                top, bottom = value[a], value[b]
                if top is not None or bottom is not None:
                    pairs.append([top, None if bottom is None else len(motifs[0])+bottom])
            rows[name] = {'name':name, 'motifs':motifs, 'cuts':pairs,
                          'suppliers':list(value['suppl']), 'rebase_id':value['id']}
    if len(rows) < 1000 or 'EcoRI' not in rows or 'Permission to use, copy, modify' not in license_text:
        raise ValueError('Unexpected source or license identity')
    return {'version':1, 'source':SOURCE, 'revision':REVISION,
            'release':'Biopython 1.88', 'rebase_release':'404 (2024)',
            'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
            'copyright':'Copyright (C) 2004. Frederic Sohm. Biopython contributors.',
            'license_name':'Biopython License Agreement', 'license_url':LICENSE,
            'license_text':license_text,
            'coordinate_convention':'Cut offsets are zero-based boundaries from the forward recognition start; each pair names top then bottom. Null in a defined pair means an uncut strand. An empty cut list means unassigned cleavage geometry.',
            'enzymes':[rows[name] for name in sorted(rows)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, help='Local original dictionary for offline conversion')
    parser.add_argument('--license', type=Path, help='Local original license for offline conversion')
    args = parser.parse_args()
    source = args.source.read_text() if args.source else urllib.request.urlopen(SOURCE, timeout=30).read().decode()
    license_text = args.license.read_text() if args.license else urllib.request.urlopen(LICENSE, timeout=30).read().decode()
    result = catalogue(source, license_text)
    output = ROOT/'catalogue'/'restriction_enzymes.json'
    output.write_text(json.dumps(result, ensure_ascii=False, separators=(',',':'))+'\n')
    print(f'{len(result["enzymes"])} enzymes → {output}')


if __name__ == '__main__':
    main()
