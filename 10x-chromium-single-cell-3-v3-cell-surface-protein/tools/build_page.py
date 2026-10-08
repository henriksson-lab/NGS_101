#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_csp as C
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-3prime-v3-cell-surface-protein.html"
def sequencing(lib):
    rows=[]
    for p in C.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Nextera Read 1 reports the 16-base cell barcode and 12-base UMI. TruSeq Read 2 reaches the 15-base antibody feature barcode after 10 diversity bases; i7 identifies the sample.')+sp.diagram(lib,C.SEQ_PRIMERS)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Cycles","Content"),C.read_layout(),scroll=False)
def render():
    product=C.feature_cdna(); lib=C.final_library()
    return "\n".join([head("10x 3' v3 Cell Surface Protein chemistry"),'<div class="wrap">','<h1>10x Chromium Single Cell 3&#39; v3 &mdash; Cell Surface Protein</h1>',info('Primary chemistry source: <a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/CG000185_ChromiumSingleCell3__FeatureBarcode_CellSurfaceProtein_Rev_B.pdf">10x Genomics User Guide CG000185 Rev B</a>.'),
      '<h2>Key oligos</h2><seq>'+oligo("Capture Sequence 1 gel-bead primer",C.bead_primer(),mods="bead")+oligo("antibody Feature Barcode oligo",C.antibody_oligo(),mods="5' antibody conjugate")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture antibody tags inside each GEM</h3>'+panel(C.capture_scene().rows(),cls="small",caption="Capture Sequence 1 pairs the antibody tag with a bead primer carrying the cell barcode and UMI; extension copies the feature barcode and Read 2 handle."),
      '<h3>(2) Co-amplify the short feature cDNA, then size-separate it from gene-expression cDNA</h3>'+panel([*Scene.duplex(list(product),label="feature cDNA").rows(),*annotation_rows(product)],cls="small",caption="Feature cDNA Primers 2 amplify between Partial Read 1N and Partial Read 2."),
      '<h3>(3) Add P5, P7 and the sample index by Feature SI PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final Cell Surface Protein library from CG000185."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
