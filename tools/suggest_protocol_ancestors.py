#!/usr/bin/env python3
"""Suggest earlier chemistry analogues for source review; never create graph edges."""
import argparse
import csv
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'lib'), str(ROOT/'catalogue/tools')]
from protocol_ancestry_similarity import rank_candidates, publication_years, cited_protocols
from properties import load


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('protocol', help='Catalogue directory ID')
    parser.add_argument('--limit', type=int, default=10)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--source', type=Path, help='Selected original text twin after manifest/version verification and mandatory primer scan; return known-DOI line pointers, not ancestry conclusions')
    args = parser.parse_args()
    rows = list(csv.DictReader((ROOT/'catalogue/ours.tsv').open(), delimiter='\t'))
    known = {row['dir']: row for row in rows}
    if args.protocol not in known:
        parser.error('Unknown protocol directory: ' + args.protocol)
    metadata = list(csv.DictReader((ROOT/'catalogue/scg_lib_structs.tsv').open(), delimiter='\t'))
    years = publication_years(rows, metadata)
    component_path = ROOT/'assets/sequence_components.json'
    components = json.loads(component_path.read_text()) if component_path.exists() else {'components': []}
    result = rank_candidates(args.protocol, rows, load(set(known)), years, components, limit=args.limit)
    if args.source:
        result['source'] = str(args.source)
        result['cited_protocols'] = cited_protocols(args.protocol, rows, args.source.read_text())
    output = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
