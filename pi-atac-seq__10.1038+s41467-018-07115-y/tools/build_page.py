#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import pi_atac as P
import seqprimers as sp
from chemdraw import Scene,annotation_rows,panel
from page import caveat,head,info,table
OUT=HERE.parent/"pi-atac-seq.html"
def main():
 lib=P.final_library(); rows=[]
 for p,h in P.primer_landings(): rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
 page="\n".join([head("Pi-ATAC-seq library chemistry"),'<div class="wrap">','<h1>Pi-ATAC-seq &mdash; joint protein indexing and chromatin accessibility</h1>',info('Defining source: <a href="https://doi.org/10.1038/s41467-018-07115-y">Chen et al., <i>Nature Communications</i> (2018)</a>.'),
 caveat('<b>Inferred adapter structure.</b> The paper reports 96 × 90 barcoding-primer combinations but does not publish their sequences. Dotted adapter regions show the standard Nextera-compatible reconstruction.'),
 '<h2>Library construction</h2><h3>(1) Fix, antibody-stain, and tagment cells in bulk</h3>'+panel(P.tagmented_scene().rows(),cls="small",caption="INFERRED — standard ATAC s5/s7 geometry. Bound Tn5 keeps the tagged DNA inside each fixed cell for sorting."),
 '<h3>(2) Index-sort one cell per well and reverse crosslinks</h3>'+info('Protein abundance is recorded by FACS; it is not encoded in the DNA library.'),
 '<h3>(3) Gap-fill and amplify with the well’s barcoding-primer pair</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="INFERRED — Nextera-compatible final library, P5 to P7'. Exact Pi-ATAC primer and index bases are unpublished."),
 '<h2>Sequencing</h2>'+info('Paired-end 2 × 75 with two per-well index reads. Primer landing sites below are computed for the inferred Nextera-compatible structure.')+sp.diagram(lib,P.SEQ_PRIMERS)+table(("Read","Primer","Primer site","First bases"),rows),'</div>'])
 OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
