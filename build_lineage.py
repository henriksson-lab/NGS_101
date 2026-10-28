#!/usr/bin/env python3
"""Build protocol-relationship data and the offline interactive graph tab."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from page import head
from protocol_lineage import build_payload


def render_page():
    return '<!doctype html><html lang="en">'+head('Protocol relationships · NGS 101',home='index.html')+'''
<link rel="stylesheet" href="assets/protocol_lineage.css">
<main class="wrap lineage-tool"><h1>Protocol relationships</h1>
<div class="lineage-controls"><label>Protocol<input id="lineage-search" type="search" list="lineage-protocols" placeholder="Find a protocol"></label><datalist id="lineage-protocols"></datalist><button id="lineage-find" type="button">Find</button>
<label><input id="lineage-documented" type="checkbox" checked> Source-documented</label><label><input id="lineage-inferred" type="checkbox" checked> Inferred</label><label><input id="lineage-similar" type="checkbox"> Sequence similarity</label><label><input id="lineage-isolated" type="checkbox"> Include unlinked</label>
<button id="lineage-focus" type="button" disabled>Focus selected</button><button id="lineage-all" type="button">All relationships</button><button id="lineage-fit" type="button">Fit</button><button id="lineage-plus" type="button" aria-label="Zoom in">+</button><button id="lineage-minus" type="button" aria-label="Zoom out">−</button><button id="lineage-export" type="button" disabled>Export SVG</button></div>
<p id="lineage-status" role="status" aria-live="polite">Loading relationships…</p>
<div class="lineage-workspace"><div id="lineage-canvas" tabindex="0" aria-label="Protocol relationship graph, drag to pan and use plus or minus to zoom"></div><aside id="lineage-detail"><p>Select a protocol or link.</p></aside></div>
<details class="lineage-method"><summary>Evidence and similarity</summary><p>Arrows show curated method development or explicitly scoped component adaptation. Solid teal: source-documented. Dashed amber: inferred. Multiple parents are allowed. Missing links mean uncurated relationships, not independent origins. A relationship never supplies another protocol’s oligos or reaction chemistry.</p>
<p>Dotted purple links are an independent sequence-similarity analysis, not ancestry. Exact known primer, adapter and fixed sequence components are deduplicated on either strand. Inferred sequences, ambiguous sequences, short components (&lt;16 nt) and low-complexity components with fewer than three base types are excluded. Similarity is IDF-weighted Jaccard: component weight ln((N+1)/(df+1)), where N is the number of protocols with eligible components and df is the component’s protocol frequency. Ubiquitous parts have zero weight. At least two components and score ≥0.15 are required; up to three nearest neighbours per protocol are retained (symmetrized). Coverage differences affect scores; no exported components means unknown coverage.</p>
<p>Only successfully built, finished, public protocols with non-illustrative declarations are shown. Source evidence and inferred rationale are available by selecting a link.</p></details>
<noscript><p>This graph viewer requires JavaScript.</p></noscript></main><script src="assets/protocol_lineage.js" defer></script></html>
'''


def build(output=ROOT,eligible=None):
    rows=list(csv.DictReader((ROOT/'catalogue/ours.tsv').open(),delimiter='\t'))
    if eligible is None:
        import build_index
        eligible=[]
        for row in rows:
            if (row.get('section') != 'published' or row.get('status') not in {'documented','construct','assay'}
                    or row.get('representation') == 'illustrative_example'):
                continue
            path=build_index.diagram_out(row['dir'])
            if path and path.exists():
                eligible.append({'id':row['dir'],'url':path.relative_to(ROOT).as_posix()})
    registry=json.loads((ROOT/'catalogue/protocol_lineage.json').read_text())
    components=json.loads((ROOT/'assets/sequence_components.json').read_text())
    payload=build_payload(rows,eligible,registry,components)
    assets=output/'assets';assets.mkdir(parents=True,exist_ok=True)
    (assets/'protocol_lineage.json').write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':'))+'\n')
    (output/'lineage.html').write_text(render_page())
    return payload


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eligible',type=Path)
    args=parser.parse_args()
    eligible=json.loads(args.eligible.read_text()) if args.eligible else None
    result=build(eligible=eligible)
    print('Protocol relationships: '+', '.join(f'{key}={value}' for key,value in result['stats'].items()))
