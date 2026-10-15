#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v3 as V
from protocol_graph import AnnotationSlot
from protocol_views import PrimerBindingSection, assemble_page
from chemdraw import Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"10x-3prime-v3.html"
def render():
    lib=V.final_library(); cdna=V.amplified_cdna(); graph=V.protocol_graph()
    return assemble_page(graph, [head("10x Chromium 3' Gene Expression v3 chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>10x Chromium Single Cell 3&#39; Gene Expression v3</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support">10x Genomics User Guide CG000183 Rev A</a>.'),
      AnnotationSlot.INTRO,
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Key oligos</h2><seq>'+oligo("gel-bead oligo",V.bead_oligo())+oligo("template-switch oligo",V.tso(),mods="3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) Capture and barcode mRNA during GEM reverse transcription</h3>'+panel(V.capture_scene().rows(),cls="small",caption="Each molecule receives a 16-base cell barcode and 12-base UMI."),
      '<h3>(2) Template-switch, break GEMs, and amplify cDNA</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The TSO supplies the second amplification handle.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Amplified cDNA retains barcode and UMI at the transcript's 3' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The source-listed adapter has a 12-bp stem and 3'-T ligation overhang."),
      '<h3>(4) Select the barcode-bearing end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final single-index v3 library. PCR restores the complete Read 2 arm."),
      PrimerBindingSection("library", "Library sequencing", "Sequencing primers and reads",
          "Paired-end, single-index sequencing. Read 1 contains the cell barcode and UMI.").render(graph),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
