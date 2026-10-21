#!/usr/bin/env python3
"""Build the offline exact-match component catalogue in isolated protocol processes.

Pass protocol directories positionally to restrict the export, or --eligible with a
JSON list of directory IDs (or objects with id/title/url). Normal builds need no
source downloads. Unavailable graphs are recorded explicitly, not silently replaced.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'lib'), str(ROOT / 'tools'), str(ROOT)]
from sequence_catalogue import combine, extract_graph


def extract_protocol(name):
    from selftest_source_oligos import graph_module_paths, load
    directory = ROOT / name
    if directory.parent != ROOT or not directory.is_dir():
        raise ValueError('protocol must be a direct directory under the repository root')
    paths = graph_module_paths(directory)
    if not paths:
        return {'status': 'missing_graph', 'occurrences': [], 'templates': []}
    if len(paths) != 1:
        raise ValueError('expected one graph declaration; found ' + str(len(paths)))
    # Modules may print diagnostic messages at import/construction time. JSON must
    # remain machine-readable while those diagnostics are captured in the report.
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        module = load(paths[0])
        sys.path.insert(0, str(paths[0].parent))
        graph = module.protocol_graph()
        result = extract_graph(graph, name)
    result['status'] = 'ok'
    if captured.getvalue().strip():
        result['diagnostic'] = captured.getvalue().strip()[-2000:]
    return result


def public_url(name):
    # Reuse the site's declaration reader for wrapper and wrapper-free protocols.
    # This discovers a declared output; it does not build or publish anything.
    from build_index import diagram_out
    output = diagram_out(name)
    return output.relative_to(ROOT).as_posix() if output else ''


def build(output, protocol_ids=None, eligible=None):
    with (ROOT / 'catalogue' / 'ours.tsv').open(newline='', encoding='utf-8') as handle:
        rows = {row['dir']: row for row in csv.DictReader(handle, delimiter='\t')}
    overrides = {}
    if eligible is not None:
        data = json.loads(Path(eligible).read_text(encoding='utf-8'))
        if isinstance(data, dict):
            data = data['protocols']
        selected = []
        for entry in data:
            name = entry['id'] if isinstance(entry, dict) else entry
            selected.append(name)
            if isinstance(entry, dict):
                overrides[name] = entry
    else:
        selected = list(protocol_ids) if protocol_ids else sorted(rows)
    unknown = set(selected) - set(rows)
    if unknown:
        raise ValueError('unknown protocols: ' + ', '.join(sorted(unknown)))
    protocols, extractions, report = [], [], []
    for name in sorted(set(selected)):
        supplied = overrides.get(name, {})
        protocols.append({'id': name, 'title': supplied.get('title', rows[name]['protocol']),
                          'url': supplied.get('url', public_url(name))})
        try:
            worker = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', name],
                                    cwd=ROOT, capture_output=True, text=True, timeout=90)
            if worker.returncode:
                raise RuntimeError(worker.stderr.strip()[-1500:] or 'worker exited ' + str(worker.returncode))
            result = json.loads(worker.stdout)
            extractions.append(result)
            item = {'protocol': name, 'status': result['status'],
                    'occurrences': len(result['occurrences']), 'templates': len(result['templates'])}
            if result.get('reason'):
                item['reason'] = result['reason']
            if result.get('unrepresented_oligos'):
                item['unrepresented_oligos'] = result['unrepresented_oligos']
            if result.get('diagnostic'):
                item['diagnostic'] = result['diagnostic']
            report.append(item)
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            report.append({'protocol': name, 'status': 'error', 'reason': str(exc)})
    result = combine(protocols, extractions, report)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(',', ':')) + '\n'
    # Keep unchanged builds from producing a needless output timestamp change.
    if not target.exists() or target.read_text(encoding='utf-8') != rendered:
        target.write_text(rendered, encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('protocols', nargs='*')
    parser.add_argument('--output', type=Path, default=ROOT / 'assets' / 'sequence_components.json')
    parser.add_argument('--eligible', type=Path)
    parser.add_argument('--worker', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        try:
            result = extract_protocol(args.worker)
        except FileNotFoundError as exc:
            result = {'status': 'missing_source', 'reason': str(exc), 'occurrences': [], 'templates': []}
        print(json.dumps(result, ensure_ascii=False))
        return 0
    result = build(args.output, args.protocols, args.eligible)
    print('Sequence catalogue: ' + ', '.join(f'{v} {k}' for k, v in result['stats'].items()))
    failures = [r for r in result['report'] if r['status'] == 'error']
    for row in failures:
        print(row['protocol'] + ': ' + row['reason'], file=sys.stderr)
    return bool(failures)


if __name__ == '__main__':
    raise SystemExit(main())
