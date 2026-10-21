#!/usr/bin/env python3
"""Build the browser-local annealing tab and bundle the shared Python model.

The normal build is offline and standard-library-only. The browser loads the pinned
Pyodide runtime on first use; model source is always served from this site.
"""
from pathlib import Path
import html
import json
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lib'))
from chemdraw import revcomp
from page import head

RUNTIME_VERSION = '314.0.7'
MODULES = ('protocol_graph', 'chemdraw', 'topology', 'annealing')


def python_bundle(*, runtime_source='cdn'):
    if runtime_source not in ('cdn', 'local'):
        raise ValueError('Unknown browser runtime source.')
    return {'runtime_version': RUNTIME_VERSION, 'runtime_source': runtime_source,
            'modules': {name: (ROOT / 'lib' / (name + '.py')).read_text(encoding='utf-8') for name in MODULES}}


def render_page(*, prototype=False):
    examples = {'linear': {'first': 'TTACGTCAGATCGTAA', 'second': 'GG' + revcomp('ACGTCAGATCGT') + 'CC'},
                'padlock': {'first': 'TTACGTCAGTGGATTGCATCGAAA', 'second': revcomp('TGCATCGA') + 'GACTAACGTAGCTAGTCA' + revcomp('ACGTCAGT')}}
    return '<!doctype html><html lang="en"><meta charset="utf-8">' + head('Anneal primers · NGS 101', home='index.html', include_header=not prototype) + '''
<link rel="stylesheet" href="assets/sequence_annotation.css"><link rel="stylesheet" href="assets/anneal.css">
<main class="wrap annotation-tool anneal-tool"><h1>Anneal primers</h1>
<form id="anneal-form" class="annotation-card"><div class="anneal-inputs">
<label for="anneal-first">Sequence 1 · 5′ → 3′<textarea id="anneal-first" rows="4" spellcheck="false" required></textarea></label>
<label for="anneal-second">Sequence 2 · 5′ → 3′<textarea id="anneal-second" rows="4" spellcheck="false" required></textarea></label></div>
<div class="anneal-controls"><label>Arrangement<select id="anneal-mode"><option value="both">Linear and padlock</option><option value="linear">Linear</option><option value="padlock">Padlock</option></select></label>
<label>Minimum paired bases (per arm)<input id="anneal-minimum" type="number" value="8" min="2" max="100" required></label>
<button id="anneal-search" type="submit" class="primary-button">Find pairings</button><button id="anneal-cancel" type="button" disabled>Cancel</button>
<button id="anneal-example" type="button" data-examples="''' + html.escape(json.dumps(examples), quote=True) + '''">Set example sequences</button></div></form>
<p id="anneal-status" role="status" aria-live="polite">Python loads on first use (usually a few seconds).</p>
<section id="anneal-results" class="annotation-card" hidden><p id="anneal-count"></p><label for="anneal-candidate">Arrangement</label><select id="anneal-candidate"></select><p id="anneal-summary"></p><div id="anneal-diagram"></div><details><summary>Coordinates</summary><pre id="anneal-coordinate-text"></pre></details><button id="anneal-download" type="button" hidden>Download SVG</button></section>
<details><summary>Scope</summary><p>Exact antiparallel A/C/G/T spans, up to 3,000 bases per input. N is allowed in unpaired regions and interrupts exact matching. Padlocks use terminal arms on the other supplied strand, trying either strand as target. Maximal arm boundaries and alternative placements remain selectable.</p><p>These are possible pairing arrangements, not thermodynamic fold predictions or evidence of binding, extension or ligation. Positions refer to the normalized 5′ → 3′ inputs. Non-selected regions in linear views are unpaired. Search limits are reported.</p></details>
<noscript><p>This tool requires JavaScript.</p></noscript></main><script src="assets/anneal_ui.js" defer></script></html>
'''



def build(output=ROOT):
    assets = output / 'assets'
    assets.mkdir(parents=True, exist_ok=True)
    (output / 'anneal.html').write_text(render_page(), encoding='utf-8')
    bundle = json.dumps(python_bundle(), ensure_ascii=False, sort_keys=True, separators=(',', ':')) + '\n'
    destination = assets / 'anneal_python.json'
    if not destination.exists() or destination.read_text(encoding='utf-8') != bundle:
        destination.write_text(bundle, encoding='utf-8')
    return output / 'anneal.html'


if __name__ == '__main__':
    print(build().relative_to(ROOT))
