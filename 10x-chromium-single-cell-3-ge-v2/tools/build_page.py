#!/usr/bin/env python3
from __future__ import annotations
import html, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v2 as V
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-3prime-v2.html"

def sequencing(lib):
    rows=[]
    for p in V.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Read 1 is 26 cycles: 16-base cell barcode then 10-base UMI. Index 1 identifies the sample; Read 2 sequences cDNA.')+table(("Read","Primer","Primer site","First bases"),rows)+table(("Read 1 cycles","Content"),V.read1_layout(),scroll=False)

def render():
    lib=V.final_library(); cdna=V.amplified_cdna()
    return "\n".join([head("10x Chromium 3' Gene Expression v2 chemistry"),'<div class="wrap">',
      '<h1>10x Chromium Single Cell 3&#39; Gene Expression v2</h1>',
      info('Primary chemistry source: <a href="https://www.10xgenomics.com/support">10x Genomics Technical Note CG000108 Rev A</a>.'),
      '<h2>Key oligos</h2><seq>'+oligo("gel-bead oligo",V.bead_oligo())+oligo("template-switch oligo",V.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture and barcode mRNA during GEM reverse transcription</h3>'+panel(V.capture_scene().rows(),cls="small",caption="Each dissolved gel bead supplies one 16-base cell barcode; each captured molecule receives a 10-base UMI."),
      '<h3>(2) Template-switch, break GEMs, and amplify cDNA</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The TSO supplies the second cDNA-amplification handle.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Amplified full-length cDNA retains barcode and UMI at the transcript's 3' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The short strand forms a 12-bp stem and leaves a 3'-T overhang for the dA-tailed cDNA fragment."),
      '<h3>(4) Select the barcode-bearing end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final v2 library. P5-side priming selects the fragment that retains Partial Read 1, cell barcode and UMI."),sequencing(lib),'</div>'])

def main():
    OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
