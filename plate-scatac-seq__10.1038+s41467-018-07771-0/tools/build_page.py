#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import plate_scatac as P
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"plate-scatac-seq.html"
def main():
 lib=P.final_library(); rows=[]
 for p,h in P.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
 page="\n".join([head("Plate_scATAC-seq library chemistry"),'<div class="wrap">','<h1>Plate_scATAC-seq &mdash; bulk tagmentation, one nucleus per PCR well</h1>',info('Defining source: <a href="https://doi.org/10.1038/s41467-018-07771-0">Chen et al., <i>Nature Communications</i> (2018)</a>.'),
 '<h2>Representative index primers</h2><seq>'+oligo("N701",P.n701())+oligo("S502",P.s502())+'</seq>',
 '<h2>Library construction</h2><h3>(1) Tagment native nuclei in bulk</h3>'+panel(P.tagmented_scene().rows(),cls="small",caption="Tn5 remains bound, holding the tagged nuclear DNA together for single-nucleus sorting."),
 '<h3>(2) Sort one nucleus into each pre-indexed well; release Tn5 and fill the gaps</h3>'+info('SDS + proteinase K releases Tn5; Tween-20 quenches SDS before the 72 °C gap-fill step.'),
 '<h3>(3) Amplify with that well’s S5xx/N7xx pair</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Plate_scATAC-seq library, P5 to P7'. Cell identity is the i5/i7 pair."),
 '<h2>Sequencing</h2>'+info('Paired-end Nextera sequencing with two 8-base index reads; the paper does not state the genomic read lengths.')+table(("Read","Primer","Primer site","First bases"),rows),'</div>'])
 OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
