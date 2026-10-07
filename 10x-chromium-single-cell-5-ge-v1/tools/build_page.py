#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v1 as V
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-5prime-v1.html"
def sequencing(lib):
    rows=[]
    for p in V.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Read 1 reports the 16-base cell barcode and 10-base UMI; i7 identifies the sample; Read 2 sequences the transcript-facing fragment.')+table(("Read","Primer","Primer site","First bases"),rows)+table(("Read 1 cycles","Content"),V.read1_layout(),scroll=False)
def render():
    lib=V.final_library(); cdna=V.amplified_cdna()
    return "\n".join([head("10x Chromium 5' Gene Expression v1 chemistry"),'<div class="wrap">','<h1>10x Chromium Single Cell 5&#39; Gene Expression v1</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support">10x Genomics Technical Note CG000109 Rev D</a>.'),
      '<h2>Key oligos</h2><seq>'+oligo("barcoded gel-bead TSO",V.bead_tso(),mods="bead; 3' rGrGrG")+oligo("poly(dT) RT primer",V.rt_primer())+'</seq>',
      '<h2>Library construction</h2><h3>(1) Prime reverse transcription from the mRNA poly(A) tail</h3>'+panel(V.capture_scene().rows(),cls="small",caption="The universal RT primer starts first-strand synthesis at the transcript 3' end."),
      '<h3>(2) Template-switch onto the barcoded gel-bead oligo</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The bead TSO writes the cell barcode and UMI at the transcript 5' end.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Bulk cDNA amplification preserves the bead-derived 5' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The short strand forms a 12-bp stem and leaves a 3'-T ligation overhang."),
      '<h3>(4) Select the bead end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final single-index v1 Gene Expression library, P5 to P7'."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
