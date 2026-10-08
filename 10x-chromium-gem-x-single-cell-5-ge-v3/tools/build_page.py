#!/usr/bin/env python3
from __future__ import annotations
import html,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v3 as V
import seqprimers as sp
from chemdraw import Scene,annotation_rows,oligo,panel
from page import head,info,table
OUT=HERE.parent/"10x-5prime-v3.html"
def sequencing(lib):
    rows=[]
    for p in V.SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} has no site")
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return '<h2>Sequencing</h2>'+info('Paired-end, dual-index sequencing: Read 1 28 cycles, i7 10 cycles, i5 10 cycles, Read 2 90 cycles. Read 1 contains the 16-base cell barcode and 12-base UMI.')+sp.diagram(lib,V.SEQ_PRIMERS)+table(("Read","Primer","Primer site","First bases"),rows)+table(("Read 1 cycles","Content"),V.read1_layout(),scroll=False)
def render():
    lib=V.final_library(); cdna=V.amplified_cdna()
    return "\n".join([head("10x GEM-X 5' Gene Expression v3 chemistry"),'<div class="wrap">','<h1>10x Chromium GEM-X Single Cell 5&#39; Gene Expression v3</h1>',info('Primary chemistry source: <a href="https://www.10xgenomics.com/support/universal-five-prime-gene-expression/documentation/steps/library-prep/chromium-gem-x-single-cell-5-v3-gene-expression-user-guide">10x Genomics User Guide CG000733 Rev A</a>.'),
      '<div class="caveat"><b>Intermediate oligos.</b> CG000733 prints the gel-bead TSO and complete final library but not the RT-primer or ligation-adapter bases. Dotted regions and steps labelled <b>INFERRED</b> reconstruct those intermediates from v1 and the final product.</div>',
      '<h2>Key oligo</h2><seq>'+oligo("barcoded GEM-X gel-bead TSO",V.bead_tso(),mods="bead; 3' rGrGrG")+'</seq>',
      '<h2>Library construction</h2><h3>(1) INFERRED &mdash; prime reverse transcription from the mRNA poly(A) tail</h3>'+panel(V.capture_scene().rows(),cls="small",caption="INFERRED — the guide names Poly-dT RT Primer B but does not print its bases."),
      '<h3>(2) INFERRED &mdash; template-switch onto the barcoded gel-bead oligo</h3>'+panel(V.switch_scene().rows(),cls="small",caption="INFERRED — the source-listed bead TSO writes the cell barcode and UMI at the transcript 5' end; dotted bases belong to the unpublished RT-primer side.")+panel([*Scene.duplex(list(cdna),label="amplified cDNA").rows(),*annotation_rows(cdna)],cls="small",caption="INFERRED — the bead-derived end is source-listed; the far cDNA-amplification handle is unpublished."),
      '<h3>(3) INFERRED &mdash; fragment, end-repair, A-tail, and ligate the Read 2 adapter</h3>'+panel(V.adapter_scene().rows(),cls="small",caption="INFERRED — adapter geometry reconstructed from the complete final product printed in CG000733."),
      '<h3>(4) Add dual sample indexes by PCR</h3>'+panel([*Scene.duplex(list(lib),label="library").rows(),*annotation_rows(lib)],cls="small",caption="Final GEM-X v3 library as printed by CG000733, with 10-base i5 and i7 sample indexes."),sequencing(lib),'</div>'])
def main(): OUT.write_text(render(),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
if __name__=="__main__": main()
