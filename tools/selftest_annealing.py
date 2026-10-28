#!/usr/bin/env python3
"""Offline behavior checks for shared exact pairing discovery and SVG export."""
from pathlib import Path
import sys
import random
import json
import shutil
import tempfile
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'lib'), str(ROOT)]
from annealing import LinearCandidate, PadlockCandidate, normalize, render, search
from chemdraw import revcomp
from checks import Check
from topology import PadlockCaptureTopology


def main():
    check = Check()
    first = 'TTACGTCAGATCGTAA'
    second = 'GG' + revcomp('ACGTCAGATCGT') + 'CC'
    result = search(first, second, 'linear')
    check('linear search returns antiparallel spans in original 5-prime coordinates',
          result['candidates'], [{'kind': 'linear', 'first_start': 2, 'second_start': 2, 'length': 12}])
    scene = LinearCandidate(2, 2, 12).scene(first, second)
    check('scene retains both full strands and their overhangs',
          [sum(len(s) for s in strand.segs) for strand in scene.strands.values()], [len(first), len(second)])
    check('second strand is drawn antiparallel by existing Scene', scene.strands['sequence 2'].rev)
    check.raises('changing a selected base cannot bypass complement validation',
                 lambda: LinearCandidate(2, 2, 12).scene(first, second[:3] + 'A' + second[4:]))
    check.raises('out-of-range span rejected', lambda: LinearCandidate(99, 2, 12).scene(first, second))
    check('FASTA and whitespace normalize without changing base meaning', normalize(' >a\nacgt \nA'), 'ACGTA')
    for bad in ('>a\nACGT\n>b\nACGT', 'ACGR', 'ACGU', 'AC-GT', '', 'A' * 3001):
        check.raises('invalid sequence rejected: ' + repr(bad[:18]), lambda bad=bad: normalize(bad))
    for kwargs in ({'minimum': True}, {'minimum': 1}, {'mode': 'automatic-chemistry'}, {'limit': 251}):
        check.raises('invalid search options rejected: ' + repr(kwargs), lambda kwargs=kwargs: search(first, second, **kwargs))
    repeated = search('ACGTCAGATCGTACGTCAGATCGT', revcomp('ACGTCAGATCGT'), 'linear')
    check('repeated exact sites remain alternatives', sorted(c['first_start'] for c in repeated['candidates']), [0, 12])
    check('no complementary sequence gives no candidates', search('A' * 20, 'A' * 20)['candidates'], [])
    check('result cap reports incomplete search', search('A' * 60, 'T' * 60, 'linear', 4, 2)['truncated'])
    target = 'TTACGTCAGTGGATTGCATCGAAA'
    probe = revcomp('TGCATCGA') + 'GACTAACGTAGCTAGTCA' + revcomp('ACGTCAGT')
    result = search(target, probe, 'padlock')
    expected = {'target': 1, 'extension_start': 2, 'extension_length': 8,
                'ligation_start': 14, 'ligation_length': 8, 'kind': 'padlock'}
    check('padlock discovers terminal arms and intervening backbone', result['candidates'], [expected])
    swapped = search(probe, target, 'padlock')
    check('either input can be the padlock', swapped['candidates'], [dict(expected, target=2)])
    model = PadlockCandidate(1, 2, 8, 14, 8).topology(target, probe)
    check('padlock uses unchanged probe arm/backbone sequence', model.ligation_arm + model.backbone + model.extension_arm, probe)
    check('bound target gap retains exact coordinates', (model.extension_site + len(model.extension_arm), model.ligation_site), (10, 14))
    repeated_target = target + 'GG' + target
    multi = search(repeated_target, probe, 'padlock')
    check('padlock repeated arm sites enumerate alternative placements', len(multi['candidates']) > 1)
    check('explicit site makes repeated padlock match renderable through shared topology',
          all(PadlockCandidate(**{k:v for k,v in c.items() if k != 'kind'}).topology(repeated_target, probe) for c in multi['candidates']))
    check.raises('implicit repeated sites remain ambiguous for existing callers',
                 lambda: PadlockCaptureTopology('ambiguous', repeated_target, model.ligation_arm, model.backbone, model.extension_arm))
    check.raises('explicit incorrect arm site is rejected',
                 lambda: PadlockCaptureTopology('bad', target, model.ligation_arm, model.backbone, model.extension_arm, ligation_site=1))
    check.raises('overlapping arm coordinates are rejected',
                 lambda: PadlockCaptureTopology('overlap', 'ACGTACGTACGT', 'ACGTACGT', 'AGTC', 'ACGTACGT', ligation_site=0, extension_site=4))
    for candidate, a, b in [(result['candidates'][0], target, probe), (search(first, second, 'linear')['candidates'][0], first, second)]:
        export = render(a, b, candidate)
        svg = ET.fromstring(export['svg'])
        check('standalone SVG has proper namespace and embedded styling', svg.tag == '{http://www.w3.org/2000/svg}svg' and 'style' in export['svg'])
        check('export remains vector and includes base text', '<text' in export['svg'] or ':text' in export['svg'])
    coordinates = render(target, probe, expected)['coordinates']
    check('render response preserves target-gap and probe-backbone coordinates',
          (coordinates['target_gap'], coordinates['probe_backbone']), ([10, 14], [8, 26]))
    check('swapped inputs keep target/probe identity and target coordinates',
          render(probe, target, dict(expected, target=2))['coordinates']['target_sequence'], 2)
    check('linear render preserves original-strand spans and overhang lengths',
          render(first, second, search(first, second, 'linear')['candidates'][0])['coordinates'],
          {'first': [2, 14], 'second': [2, 14], 'first_overhangs': [2, 2], 'second_overhangs': [2, 2]})
    unknown_probe = model.ligation_arm + 'NNNN' + model.extension_arm
    unknown = search(target, unknown_probe, 'padlock')
    check('unknown backbone bases are retained without preventing exact terminal arms', len(unknown['candidates']), 1)
    unknown_linear = search('AAAAANNNNN', 'NNNNNTTTTT', 'linear', 2)['candidates']
    check('unknown bases are never counted as exact linear pairs',
          max(c['length'] for c in unknown_linear) == 5 and all(c['first_start'] + c['length'] <= 5 and c['second_start'] >= 5 for c in unknown_linear))
    check('all-unknown inputs cannot create a pairing', search('N' * 20, 'N' * 20)['candidates'], [])
    check.raises('unknown selected span cannot bypass exact-pair check', lambda: LinearCandidate(0, 0, 4).scene('NNNN', 'NNNN'))
    check('unknown padlock backbone renders with the shared topology', '4 backbone bases' in render(target, unknown_probe, unknown['candidates'][0])['summary'])
    # Independent short-sequence oracle enumerates arm lengths before filtering
    # dominated boundaries, rather than deriving lengths from terminal seed matches.
    rng = random.Random(27)
    oracle_ok = True
    for _ in range(12):
        a = ''.join(rng.choice('ACGT') for _ in range(10))
        b = revcomp(a[6:8]) + 'GA' + revcomp(a[1:3])
        expected_all = set()
        for number, target_seq, probe_seq in ((1, a, b), (2, b, a)):
            raw = []
            for e in range(2, len(probe_seq) - 2):
                for l in range(2, len(probe_seq) - e):
                    for es in range(len(target_seq) - e + 1):
                        if target_seq[es:es + e] != revcomp(probe_seq[-e:]):
                            continue
                        for ls in range(len(target_seq) - l + 1):
                            if target_seq[ls:ls + l] == revcomp(probe_seq[:l]) and (es + e <= ls or ls + l <= es):
                                raw.append((es, e, ls, l))
            for es, e, ls, l in raw:
                if not any(xs == es and ys + yl == ls + l and xe >= e and yl >= l and (xe > e or yl > l)
                           for xs, xe, ys, yl in raw):
                    expected_all.add((number, es, e, ls, l))
        observed = search(a, b, 'padlock', 2)
        got = {(c['target'], c['extension_start'], c['extension_length'], c['ligation_start'], c['ligation_length']) for c in observed['candidates']}
        oracle_ok = oracle_ok and not observed['truncated'] and got == expected_all
    check('maximal padlock discovery agrees with exhaustive arm-length oracle', oracle_ok)
    check.raises('unknown candidate type is rejected', lambda: render(first, second, {'kind': 'thermodynamic'}))
    # Exercise the normal builder/assembler, including dependencies that HTML
    # href/src discovery cannot see (Worker() and fetch() targets).
    import build_anneal
    import build_site
    with tempfile.TemporaryDirectory(prefix='chem-anneal-build-test-') as temp:
        root = Path(temp)
        build_anneal.build(root)
        bundle = json.loads((root / 'assets' / 'anneal_python.json').read_text())
        check('production bundle uses the pinned CDN configuration without fetching runtime files',
              bundle['runtime_source'] == 'cdn' and not (root / 'assets' / 'pyodide').exists())
        check('bundle contains the fresh shared Python sources rather than another implementation',
              set(bundle['modules']) == set(build_anneal.MODULES) and
              all(source == (ROOT / 'lib' / (name + '.py')).read_text() for name, source in bundle['modules'].items()))
        page = (root / 'anneal.html').read_text()
        check('production page includes shared navigation and no prototype restriction',
              'href="anneal.html"' in page and 'prototype' not in page.lower() and 'pending' not in page.lower())
        for name in ('index', 'about', 'annotate', 'samplesheet', 'restriction', 'tm', 'source-downloads'):
            (root / (name + '.html')).write_text('<!doctype html><title>Build fixture</title>')
        (root / 'assets' / 'sequence_components.json').write_text('{}')
        (root / 'assets' / 'restriction_catalogue.json').write_text('{}')
        (root / 'assets' / 'restriction_worker.js').write_text('')
        for name in ('anneal_worker.mjs', 'anneal_ui.js', 'anneal.css', 'sequence_annotation.css'):
            shutil.copyfile(ROOT / 'assets' / name, root / 'assets' / name)
        original_root, original_site = build_site.ROOT, build_site.SITE
        try:
            build_site.ROOT, build_site.SITE = root, root / '_site'
            build_site.assemble({})
            check('site assembly includes all annealing worker, JSON, UI and CSS dependencies',
                  all((build_site.SITE / path).is_file() for path in (
                      'anneal.html', 'assets/anneal_worker.mjs', 'assets/anneal_python.json',
                      'assets/anneal_ui.js', 'assets/anneal.css', 'assets/sequence_annotation.css')))
            check('assembled annealing navigation and stylesheet/script links resolve', build_site.check_links(), [])
        finally:
            build_site.ROOT, build_site.SITE = original_root, original_site
    return check.report()


if __name__ == '__main__':
    main()
