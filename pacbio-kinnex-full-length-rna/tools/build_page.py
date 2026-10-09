#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import kinnex as K
from chemdraw import Row, annotation_rows, panel, strand_row
from dumbbell import dumbbell_rows
from page import head, info
OUT=H.parent/"kinnex-full-length-rna.html"

def render():
    cycle=K.polymerase_cycle()
    parts=[head(K.TITLE),'<div class="wrap">',f'<h1>{K.TITLE}</h1>',
           '<p class="research-notes"><a href="01_kinnex-full-length-rna.html">Research notes</a></p>',
           info(K.SOURCE),info(K.SUMMARY),f'<div class="caveat">{K.CAVEAT}</div>']
    for h,rows,cap in K.sections(): parts.extend([f'<h2>{h}</h2>',panel(rows,cls="long",caption=cap)])
    parts.extend(['<h2>Sequencing-primer and polymerase binding</h2>',
        panel([Row(chunks=[("sequencing primer + polymerase  --->  [terminal-adapter primer site]", "r2", False)]),
               *dumbbell_rows(K.SMRTBELL)],cls="long",caption="The sequencing primer anneals in a proprietary terminal adapter; its bases are unavailable, so the binding site and extension direction are shown without a fabricated sequence."),
        '<h2>HiFi passes and segmentation</h2>',
        panel([strand_row(cycle),*annotation_rows(cycle)],cls="long",caption="One polymerase circuit traverses the eight-insert array in both orientations. Consensus is formed first; the seven ordered segmentation junctions then delimit the eight original cDNAs."),'</div>'])
    return '\n'.join(parts)

if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
