#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v2 as V
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"10x-3prime-v2.html"

def sequencing(lib):
    return sp.section(lib,V.SEQ_PRIMERS,heading="Sequencing",intro='Read 1 is 26 cycles: 16-base cell barcode then 10-base UMI. Index 1 identifies the sample; Read 2 sequences cDNA.',required_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths={"Read 1":26,"Index 1 (i7)":8},show_read_content=True)

def render():
    lib=V.final_library(); cdna=V.amplified_cdna(); graph=V.protocol_graph()
    return "\n".join([head("10x Chromium 3' Gene Expression v2 chemistry"),'<div class="wrap" data-protocol-complete="true">',
      '<h1>10x Chromium Single Cell 3&#39; Gene Expression v2</h1>',
      info('Primary chemistry source: <a href="https://www.10xgenomics.com/support">10x Genomics Technical Note CG000108 Rev A</a>.'),
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Key oligos</h2><seq>'+oligo("gel-bead oligo",V.bead_oligo())+oligo("template-switch oligo",V.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture and barcode mRNA during GEM reverse transcription</h3>'+panel(V.capture_scene().rows(),cls="small",caption="Each dissolved gel bead supplies one 16-base cell barcode; each captured molecule receives a 10-base UMI."),
      '<h3>(2) Template-switch, break GEMs, and amplify cDNA</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The TSO supplies the second cDNA-amplification handle.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Amplified full-length cDNA retains barcode and UMI at the transcript's 3' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The short strand forms a 12-bp stem and leaves a 3'-T overhang for the dA-tailed cDNA fragment."),
      '<h3>(4) Select the barcode-bearing end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final v2 library. P5-side priming selects the fragment that retains Partial Read 1, cell barcode and UMI."),sequencing(lib),'</div>'])

def main():
    OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
