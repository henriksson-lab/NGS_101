#!/usr/bin/env python3
"""Build an opt-in annealing review prototype outside the checkout.

Requires a separately downloaded Pyodide core directory. This standalone benchmark builder never changes normal
site output and does not fetch or install third-party dependencies.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
sys.path.insert(0, str(ROOT))
from build_anneal import RUNTIME_VERSION, python_bundle, render_page

RUNTIME_FILES = ('pyodide.mjs', 'pyodide.asm.mjs', 'pyodide.asm.wasm', 'python_stdlib.zip', 'pyodide-lock.json')


def build(output: Path, runtime: Path, cdn_runtime: bool = False):
    output = output.resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError('Use build_anneal.py for site output; keep the standalone benchmark outside the checkout.')
    missing = [name for name in RUNTIME_FILES if not (runtime / name).is_file()]
    if missing:
        raise ValueError('Incomplete Pyodide core: ' + ', '.join(missing))
    assets = output / 'assets'
    target = assets / 'pyodide'
    target.mkdir(parents=True, exist_ok=True)
    files = []
    for name in RUNTIME_FILES:
        source = runtime / name
        if not cdn_runtime:
            shutil.copyfile(source, target / name)
        files.append({'name': name, 'bytes': source.stat().st_size, 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()})
    for name in ('anneal_worker.mjs', 'anneal_ui.js', 'anneal.css', 'sequence_annotation.css'):
        shutil.copyfile(ROOT / 'assets' / name, assets / name)
    (assets / 'anneal_python.json').write_text(json.dumps(python_bundle(runtime_source='cdn' if cdn_runtime else 'local')))
    (output / 'anneal.html').write_text(render_page(prototype=True))
    (output / 'runtime-manifest.json').write_text(json.dumps({'version': RUNTIME_VERSION, 'source': 'https://cdn.jsdelivr.net/pyodide/v' + RUNTIME_VERSION + '/full/', 'files': files}, indent=2))
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--cdn-runtime', action='store_true', help='Load the pinned runtime from its CDN in the opt-in prototype; model files stay local.')
    args = parser.parse_args()
    print(build(args.output, args.runtime, args.cdn_runtime))
