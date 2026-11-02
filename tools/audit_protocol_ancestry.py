#!/usr/bin/env python3
"""List missing predecessor links and ancestry review outcomes (offline)."""
import argparse
import csv
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
from protocol_ancestry import ancestry_audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--filter', choices=['no-parent', 'no-formal-parent', 'unreviewed', 'blocked', 'external', 'cannot-assign', 'inferred', 'guessed', 'unlinked', 'all'], default='no-parent')
    parser.add_argument('--public', action='store_true', help='Only protocols in the generated public graph; still count hidden registry parents')
    parser.add_argument('--format', choices=['tsv', 'json'], default='tsv')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--summary', action='store_true', help='Print aggregate counts instead of protocol rows')
    parser.add_argument('--check-reviewed', action='store_true', help='Fail if any selected-scope protocol is still unreviewed; explicit source blockers count as reviewed outcomes')
    args = parser.parse_args()
    rows = list(csv.DictReader((ROOT / 'catalogue/ours.tsv').open(), delimiter='\t'))
    registry = json.loads((ROOT / 'catalogue/protocol_lineage.json').read_text())
    reviews = json.loads((ROOT / 'catalogue/protocol_ancestry_reviews.json').read_text())
    public_path = ROOT / 'assets/protocol_lineage.json'
    visible = {node['id'] for node in json.loads(public_path.read_text())['nodes']} if public_path.exists() else None
    if args.public and visible is None:
        parser.error('Build the public relationship graph first: python3 build_lineage.py')
    report = ancestry_audit(rows, registry, reviews, visible_ids=visible)
    scope = [item for item in report['protocols'] if not args.public or item['id'] in visible]
    predicates = {'no-parent': lambda x: x['without_internal_parent'],
                  'no-formal-parent': lambda x: x['without_formal_parent'],
                  'unreviewed': lambda x: x['review_status'] == 'unreviewed',
                  'blocked': lambda x: x['has_source_blocker'],
                  'external': lambda x: bool(x['external_predecessors']),
                  'cannot-assign': lambda x: x['formal_assignment'] == 'cannot_formally_assign',
                  'inferred': lambda x: bool(x['inferred_parents']) or any(p['basis'] == 'inferred' for p in x['external_predecessors']),
                  'guessed': lambda x: bool(x['guessed_parents']) or any(p['basis'] == 'guessed' for p in x['external_predecessors']),
                  'unlinked': lambda x: x['unlinked'], 'all': lambda x: True}
    selected = [item for item in scope if predicates[args.filter](item)]
    if args.summary:
        from collections import Counter
        output = json.dumps({'protocols': len(scope), 'matching_filter': len(selected),
                             'without_internal_parent': sum(x['without_internal_parent'] for x in scope),
                             'without_formal_parent': sum(x['without_formal_parent'] for x in scope),
                             'review_states': dict(Counter(x['review_status'] for x in scope)),
                             'formal_assignments': dict(Counter(x['formal_assignment'] for x in scope))}, indent=2) + '\n'
    elif args.format == 'json':
        output = json.dumps({'scope': report['scope'], 'filter': args.filter, 'protocols': selected}, indent=2, ensure_ascii=False) + '\n'
    else:
        import io
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
        writer.writerow(('protocol', 'name', 'review', 'formal_assignment', 'internal_parents', 'external_predecessors', 'inference', 'reason', 'blocked_sources'))
        for item in selected:
            writer.writerow((item['id'], item['title'], item['review_status'], item['formal_assignment'], ','.join(item['internal_parents']), '; '.join(p['title'] for p in item['external_predecessors']), item['inference_status'], item['reason'], '; '.join(blocker['url'] for blocker in item['blockers'])))
        output = buffer.getvalue()
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')
    return int(args.check_reviewed and any(item['review_status'] == 'unreviewed' for item in scope))


if __name__ == '__main__':
    raise SystemExit(main())
