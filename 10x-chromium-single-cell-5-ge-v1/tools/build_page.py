#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v1 as V
from protocol_views import PrimerBindingSection,render_annotations
from protocol_graph import AnnotationSlot
from chemdraw import Segment,Scene,annotation_rows,oligo,panel,protocol_graph_panels
from page import head,info
OUT=HERE.parent/"10x-5prime-v1.html"
def render():
    lib=V.final_library(); cdna=V.amplified_cdna(); graph=V.protocol_graph()
    return "\n".join([head("10x Chromium 5' Gene Expression v1 chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>10x Chromium Single Cell 5&#39; Gene Expression v1</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support">10x Genomics Technical Note CG000109 Rev D</a> (Table 2B, 5&#39; Gene Expression). It prints oligos and products but no reaction conditions; the revision-matched user guide CG000086 it cites was not located, so enzymes and conditions are left unnamed. INFERRED: the pre-switch first strand with its untemplated CCC tail is drawn from the template-switching mechanism; CG000109 prints only the switched reverse-transcript product. Single-stranded products are drawn as duplexes for column alignment.'),
      render_annotations(graph,AnnotationSlot.INTRO),
      '<h2>Reaction workflow</h2>'+protocol_graph_panels(graph,cls="long"),
      '<h2>Key oligos</h2><seq>'+oligo("barcoded gel-bead TSO",V.bead_tso(),mods="bead; 3' rGrGrG")+oligo("Poly-dT RT Primer PN-2000007",V.rt_primer())+oligo("cDNA Primer Mix PN-220106 forward",[Segment("Partial Read 1",V.CDNA_FORWARD,"r1")])+oligo("cDNA Primer Mix PN-220106 reverse (Table 2B)",[Segment("Non-Poly(dT)",V.CDNA_REVERSE,"tso")])+oligo("Adaptor (Read 2) PN-220026 top",V.adapter_top())+oligo("Adaptor (Read 2) PN-220026 short strand",V.adapter_bottom())+oligo("Sample Index PCR Primer PN-220111",[Segment("P5 - Partial Read 1",V.INDEX_FORWARD,"p5")])+oligo("Sample Index PCR Primer PN-220103",[Segment("P7",V.INDEX_REVERSE[:24],"p7"),Segment("sample index","N"*8,"cbc",placeholder=True),Segment("Partial Read 2",V.INDEX_REVERSE[32:],"r2")])+'</seq>',
      '<h2>Library construction</h2><h3>(1) Prime reverse transcription from the mRNA poly(A) tail</h3>'+panel(V.capture_scene().rows(),cls="small",caption="The universal RT primer starts first-strand synthesis at the transcript 3' end."),
      '<h3>(2) Template-switch onto the barcoded gel-bead oligo</h3>'+panel(V.switch_scene().rows(),cls="small",caption="The bead TSO writes the cell barcode and UMI at the transcript 5' end.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="Bulk cDNA amplification preserves the bead-derived 5' end."),
      '<h3>(3) Fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="The short strand forms a 12-bp stem and leaves a 3'-T ligation overhang."),
      '<h3>(4) Select the bead end by sample-index PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final single-index v1 Gene Expression library, P5 to P7'."),
      PrimerBindingSection("library", "Library sequencing", "Sequencing",
          'CG000109 Rev D Figure 2 recommends a 26-cycle Read 1 (16-base cell barcode and 10-base UMI), 8-cycle i7, and 98-cycle Read 2 for this v1 5-prime Gene Expression library. Read 2 sequences the transcript-facing fragment.').render(graph),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
