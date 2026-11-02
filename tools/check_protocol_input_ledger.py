#!/usr/bin/env python3
"""Check finite input-review coverage against both inventoried candidate sets.

This checks review accounting and evidence locators, never the truth of a source
conclusion or biochemical compatibility. Source-unavailable rows may be closed as
blocked reviews; they do not become confirmed missing-information claims.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'catalogue/curation_rounds/2026-11-02-full-protocol-input-ledger.json'


def key(row):
    return row['scope'], row['protocol'], row['material']


def validate(ledger, public, notes):
    expected = set()
    for scope, inventory in [('public', public), ('notes', notes)]:
        for protocol in inventory['protocols']:
            for material in protocol['inputs']:
                expected.add((scope, protocol['protocol'], material['material']))
    issues = []
    rows = ledger.get('rows', [])
    keys = [key(row) for row in rows]
    if len(keys) != len(set(keys)):
        issues.append('duplicate input-review keys')
    for missing in sorted(expected - set(keys)):
        issues.append('candidate lacks review row: ' + repr(missing))
    for stale in sorted(set(keys) - expected):
        issues.append('review row not in current inventory: ' + repr(stale))
    manual = ledger.get('manual_entries', [])
    expected_manual = {p['protocol'] for p in public['protocols'] if p['status'] == 'no_unique_graph'}
    manual_ids = [row['protocol'] for row in manual]
    if len(manual_ids) != len(set(manual_ids)) or set(manual_ids) != expected_manual:
        issues.append('manual assay reviews do not match entries without a unique graph')
    for row in [*rows, *manual]:
        review = row.get('review', {})
        identifier = row['protocol'] + '/' + row['material']
        if review.get('status') not in {'reviewed', 'blocked', 'deferred'}:
            issues.append('unfinished review: ' + identifier)
        if not review.get('disposition') or not review.get('rationale') or not review.get('limits'):
            issues.append('review needs disposition, rationale and limits: ' + identifier)
        originals = review.get('originals', [])
        if review.get('status') == 'reviewed' and not originals:
            issues.append('reviewed classification needs inspected-source locators: ' + identifier)
        for original in originals:
            if not all(original.get(field) for field in ('file', 'identity', 'sections', 'availability')):
                issues.append('source needs identity, locator and availability: ' + identifier)
        if review.get('status') in {'blocked', 'deferred'} and not review.get('deferred'):
            issues.append('blocked/deferred review needs an actual target or model blocker: ' + identifier)
        if review.get('absence_claim'):
            issues.append('input review cannot accept a missing-information claim: ' + identifier)
    return issues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ledger', nargs='?', type=Path, default=DEFAULT)
    args = parser.parse_args()
    ledger = json.loads(args.ledger.read_text())
    public = json.loads((ROOT / 'catalogue/protocol_input_inventory.json').read_text())
    notes = json.loads((ROOT / 'catalogue/protocol_input_notes_candidates.json').read_text())
    problems = validate(ledger, public, notes)
    for problem in problems:
        print('FAIL', problem)
    if problems:
        return 1
    print(f"Input ledger: {len(ledger['rows'])} candidates and {len(ledger['manual_entries'])} manual assay reviews accounted for; source and model deferrals remain explicit.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
